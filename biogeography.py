from typing import Dict, Any, List
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3D projection)
import matplotlib.cm as cm
import matplotlib.image as mpimg

import paleobiodb_interface as pbdb
from paleobiodb_interface import rv
from cladogram import cladogram_lengths

def plot_taxon_occurrences_3d(
    taxon_identifier: Any,
    *,
    coordinate_filters: Dict[str, Any] = None,
    plotting_preferences: Dict[str, Any] = None,
    time_range: List[float] = None,
) -> None:
    """
    Parameters:
    - taxon_identifier: name string (e.g. "Trilobita") or PBDB taxon_no (int).
    - coordinate_filters: dict controlling coord quality thresholds.
    - plotting_preferences: colormap, marker size, alpha, figsize, dpi.
    - time_range: [min_ma, max_ma] to filter occurrences by age (Ma).
    """

    # Fetch occurrences (if internet allowed) or expect the caller to have pre-fetched data
    occurrences: List[Dict[str, Any]] = pbdb.query_occs_by_taxon(taxon_identifier)
    if not occurrences:
        raise ValueError("No occurrences returned from PBDB for the specified taxon")
        
    # Prepare taxon identifiers for distance computation
    # cladogram_distances, max_dist = cladogram_lengths(taxon_identifier)

    # Extract required fields and filter coordinates
    lats, lons, ages, cladd = ([] for i in range(4))

    maxage= 0.0
    max_dist = 0
    for occ in occurrences:
        # Skip records with missing required fields
        if rv.LAT not in occ or rv.LON not in occ or rv.MIN_MA not in occ or rv.MAX_MA not in occ or rv.TAXON_ID not in occ or rv.GENUS not in occ:
            continue

        max_ma = occ[rv.MAX_MA]
        min_ma = occ[rv.MIN_MA]
        age = (min_ma + max_ma)/2
        if age > maxage:
            maxage = age
                
        # Filter results outside of the requested timespan TODO: don't query outside the requested range
        if time_range is not None:
            min_allowed, max_allowed = float(time_range[0]), float(time_range[1])
            if age > max_allowed or age < min_allowed:
                continue
                
        # if occ[rv.TAXON_ID] not in cladogram_distances.keys():
        #     continue


        lats.append(float(occ[rv.LAT]))
        lons.append(float(occ[rv.LON]))
        ages.append(age)
        # cladd.append(cladogram_distances[occ[rv.TAXON_ID]])
        
        newhash = hash(occ[rv.GENUS]) % 100 # Dummy distance based on hash of genus name

        cladd.append(newhash)  
        if newhash > max_dist:
            max_dist = newhash
                 
    # Build arrays for plotting
    lats = np.array(lats)
    lngs = np.array(lons)
    zvals = -np.array(ages) # inverted so recent (small Ma) are higher numbers
    distances = np.array(cladd)/max_dist if max_dist > 0 else np.zeros_like(np.array(cladd))
    genera = set(cladd)

    # Prepare plotting preferences
    plotting_preferences = plotting_preferences or {}
    cmap = plotting_preferences.get('cmap', 'gist_rainbow')
    marker = plotting_preferences.get('marker', 'o')
    ms = plotting_preferences.get('markersize', 8)
    alpha = plotting_preferences.get('alpha', 0.8)
    figsize = plotting_preferences.get('figsize', (12, 6))
    dpi = plotting_preferences.get('dpi', 150)

    img = mpimg.imread('plate_carree.png')

    # Build the 3D scatter
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(lngs, lats, zvals, c=distances, cmap=cmap, s=ms, marker=marker)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Age (Ma)')
    ax.set_title(f'Occurrences for {taxon_identifier} (n={len(lats)} of {len(occurrences)}, {len(genera)} genera)')

    X1,Y1 = np.meshgrid(np.linspace(-180, 180, img.shape[1]), np.linspace(-90, 90, img.shape[0]))

    ax.plot_surface(X1, -Y1, np.zeros_like(X1) - maxage*1.05, rstride=29, cstride=29, facecolors=img, shade=False) # Works but is slow/looks bad
    ax.set_box_aspect((1, 1, 0.5)) 
    ax.set_aspect('equalxy')
    # ax.imshow(img, extent=[-180, 180, -90, 90], aspect='auto', origin='lower', interpolation='nearest') # Doesn't work
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

    # Colorbar
    # cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    # cbar.set_label('Cladogram distance')

    # Finalize: save or show
    plt.show()

def save_and_plot_platecarree_image(output_png: str = 'plate_carree.png', *, dpi: int = 150, figsize=(12, 6), show: bool = True):
    """Create a global PlateCarree map, save it as a PNG, then load and plot it

    The loaded image is plotted using geographic axis scales so x ranges from
    -180..180 (longitude) and y ranges from -90..90 (latitude).

    Parameters:
    - output_png: path to save the generated map image.
    - dpi: image resolution when saving.
    - figsize: tuple passed to matplotlib figure size.
    - show: if True, call plt.show() on the final figure.

    Returns:
    - fig, ax: the matplotlib Figure and Axes with the image plotted (ax is not a Cartopy axes).
    """

    # 1) Create and save a cartopy PlateCarree figure covering the globe
    fig1 = plt.figure(figsize=figsize, dpi=dpi)
    ax1 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    # Ensure the map covers the whole globe
    ax1.set_global()
    # ax1.coastlines(resolution='110m')
    ax1.add_feature(cfeature.LAND, facecolor='lightgray')
    ax1.add_feature(cfeature.OCEAN, facecolor='white')
    # Force the projection extent to geographic coordinates
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    # Remove whitespace/margins so the saved image maps exactly to geographic extent
    fig1.subplots_adjust(left=0, right=1, bottom=0, top=1)
    # Save to file with no padding
    fig1.savefig(output_png, bbox_inches='tight', pad_inches=0, dpi=dpi)
    plt.close(fig1)

    # 2) Load the saved image and plot it into a regular matplotlib axes

    img = mpimg.imread(output_png)

    fig2, ax2 = plt.subplots(figsize=figsize)
    # Plot the image and map its pixel extents to geographic coordinates
    ax2.imshow(img, extent=[-180, 180, -90, 90], aspect='auto', zorder=1)
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(-90, 90)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.set_title(f'Loaded map image: {output_png}')

    # Optional: overlay gridlines for verification
    ax2.set_xticks(np.linspace(-180, 180, 13))
    ax2.set_yticks(np.linspace(-90, 90, 7))
    ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)

    if show:
        plt.show()

    return fig2, ax2


def command_line_interface():
    import argparse

    parser = argparse.ArgumentParser(description='Plot occurrences of a taxon from PaleobioDB in 3D.')
    parser.add_argument('taxon', type=str, help='Taxon name or PBDB taxon_no to plot occurrences for.')
    args = parser.parse_args()
    plot_taxon_occurrences_3d(args.taxon)

if __name__ == "__main__":
    # Example usage
    command_line_interface()
    # save_and_plot_platecarree_image('plate_carree.png', dpi=150, figsize=(12,6), show=True)