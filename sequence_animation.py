from mpl_toolkits.basemap import Basemap
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt
import requests as req
import pickle
from tqdm import tqdm
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.image import BboxImage
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.legend_handler import HandlerPolyCollection
from matplotlib.collections import PolyCollection
import operator


baseUrl = 'https://macrostrat.org/api/'
stagesQuery = f'{baseUrl}columns?project_id=1'
paleoflowQuery = f'{baseUrl}measurements?measurement_id=80&project_id=1&response=long&show_values'

fname = 'column_locs.pkl'
animation_fname = 'sequence_animation_2ma.pkl'
paleoflow_fname = 'paleoflow_brandchadwick.pkl'
flow_animation_fname = 'paleoflow_animation.pkl'
out_image_fname = 'na-macrostrat.gif'
frames = 270
frame_delay = 30 #ms

plot_columns = True
plot_paleoflows = True
animate = True
save_image = False
show_plot = True

max_age = 540 # Ma
step = max_age/frames # Ma


print('Getting column location information...')
try:
    with open(fname, 'rb') as f:
        coldata = pickle.load(f)

except Exception as e:
    res = req.get(stagesQuery)
    data = res.json()['success']['data']

    coldata = {}
    for col in data:
        coldata[col['col_id']] = (float(col['lat']), float(col['lng']))
    with open(fname, 'wb') as f:
        pickle.dump(coldata, f)

print('Getting Macrostrat gap bound package data...')

try:
    with open(animation_fname, 'rb') as f:
        animation_data = pickle.load(f)

except:
    x = np.arange(0, max_age, step)
    animation_data = []
    for i, interval in enumerate(tqdm(x)):
        res = req.get(baseUrl+f'sections?age_top={interval}&age_bottom={interval+step}')
        column_ids = [x['col_id'] for x in res.json()['success']['data']]
        animation_data.append(column_ids)
    
    animation_data = animation_data[::-1] # Reverse to go in time order
    with open(animation_fname, 'wb') as f:
        pickle.dump(animation_data, f)

print('Getting Brand and Chadwick (2015) paleocurrent data...')

try:
    with open(paleoflow_fname, 'rb') as f:
        paleoflowData = pickle.load(f)

except:
    res = req.get(paleoflowQuery)
    flow_data = res.json()['success']['data']
    paleoflowData = [ dict(azimuth=x['measure_value'][0], lat=x['lat'], lon=x['lng'], err=x['measure_error'][0], unit_id=x['unit_id']) for x in flow_data]

    print('Enriching paleocurrent data with dates...')
    units_to_query = set()
    for dict in paleoflowData:
        units_to_query.add(dict['unit_id'])

    list_to_query = list(units_to_query)
    list_to_query.sort() # Sort by unit_id
    ages = []

    for unit in tqdm(list_to_query):
        res = req.get(baseUrl+f'units?unit_id={unit}')
        unit_data = res.json()['success']['data'][0]
        ages.append((unit_data['t_age'] + unit_data['b_age'])/2)

    paleoflowData.sort(key=operator.itemgetter('unit_id')) # Sort by unit_id for optimal matching

    flowPointer = 0
    for unit, age in zip(list_to_query, ages):
        while flowPointer < len(paleoflowData) and paleoflowData[flowPointer]['unit_id'] == unit:
            paleoflowData[flowPointer]['age'] = age
            flowPointer += 1

    paleoflowData.sort(key=operator.itemgetter('age')) # Sort by age for animating

    with open(paleoflow_fname, 'wb') as f:
        pickle.dump(paleoflowData, f)
    pass

print('Age binning paleoflows...')

try:
    with open(flow_animation_fname, 'rb') as f:
        flow_animation = pickle.load(f)

except:
    x = np.arange(0, max_age, step)
    flow_animation = []
    for i, interval in enumerate(x):
        flow_list = [1 if d['age'] > interval and d['age'] <= interval+step else 0 for d in paleoflowData]
        flow_animation.append(flow_list)
    
    with open(flow_animation_fname, 'wb') as f:
        pickle.dump(flow_animation, f)

flow_animation = flow_animation[::-1] # Reverse to go in time order


def extract_coords(dict):
    longs = []
    lats = []
    for lat, lng in dict.values():
        longs.append(lng)
        lats.append(lat)
    return lats, longs

def extract_arrows(list):
    longs = []
    lats = []
    for dict in list:
        longs.append(dict['lon'])
        lats.append(dict['lat'])
    return lats, longs

def extract_uv(list):
    us = []
    vs = []
    for dict in list:
        az = (90-dict['azimuth']) * 2 * np.pi/360
        us.append(np.cos(az))
        vs.append(np.sin(az))
    return us, vs

print('Plotting data...')
lats, longs = extract_coords(coldata)

plt.figure(figsize=(10.5,12))

# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS84 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(width=7000000,height=8000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=3000.,projection='lcc',\
            lat_1=35.,lat_2=55,lat_0=50,lon_0=-103.)
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.fillcontinents(color='ivory',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,10.))
m.drawmeridians(np.arange(-180.,181.,10.))
m.drawmapboundary(fill_color='aqua') 

if plot_columns:
    x, y = m(longs, lats)
    column_dots = m.scatter(x,y,50,marker='o',color='k', label='Gap bound packages')

# Custom colormap for paleoflow directions
cmap = mpl.colormaps['Set2']
# Normalization: values from 0 to 2pi
norm = Normalize(vmin=0, vmax=2*np.pi)

if plot_paleoflows:
    arrow_lat, arrow_lon = extract_arrows(paleoflowData)
    us, vs = extract_uv(paleoflowData)
    arr_x, arr_y = m(arrow_lon, arrow_lat)

    flows = plt.quiver(arr_x, arr_y, us, vs, np.arctan2(us, vs)+np.pi, cmap=cmap, norm=norm, pivot='tail', angles='xy', scale=25, label='Paleocurrents')

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='black')

    # place a text box in upper left in axes coords
    megasequence_text = plt.text(0.04, 0.17, "Macrostratigraphy", transform=plt.gca().transAxes, color='white', fontsize=28,
            verticalalignment='top', fontweight='bold', bbox=props)


plt.title(f"North American Macrostratigraphy By Time ({step} Ma bins)")

viridis = mpl.colormaps['viridis']

def update(frame):
    res = []
    if plot_columns:
        filtered_dict = {k:v for k,v in coldata.items() if k in animation_data[frame]}
        lats, longs = extract_coords(filtered_dict)
        x, y = m(longs, lats)
        data = np.stack([x, y]).T
        column_dots.set_offsets(data)
        column_dots.set_color(viridis(frame/frames))
        res.append(column_dots)

    if plot_paleoflows:
        # Indicator list: 0 for transparent, 1 for colored
        color_indicator = flow_animation[frame]
        flows.set_alpha(color_indicator)
        res.append(flows)

    age = max_age - step * frame
    megasequences = {max_age: ['Sauk', 'xkcd:snot green'], 462: ['Tippecanoe', 'xkcd:slate green'], 396: ['Kaskaskia', 'xkcd:peach'], 
                     324: ['Appalachian (Absaroka)', 'xkcd:warm purple'], 252: ['Triassic (Absaroka)', 'xkcd:pinky purple'], 186: ['Jurassic (Absaroka)', 'xkcd:aquamarine'], 
                     132: ['Zuni', 'xkcd:greyish brown'], 60: ['Tejas', 'xkcd:sunny yellow'], 24: ['post-Tejas', 'xkcd:buff']}
    boundaries = np.array([max_age, 462, 396, 324, 252, 186, 132, 60, 24])
    index = np.count_nonzero(boundaries >= age) - 1
    megasequence = megasequences[boundaries[index]]

    megasequence_text.set_text(megasequence[0])
    megasequence_text.set_color(megasequence[-1])
    res.append(megasequence_text)

    return tuple(res)

class HandlerImage(HandlerPolyCollection):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        
        x = np.linspace(1, -1, 100)
        y = np.linspace(-1, 1, 100)
        # full coordinate arrays
        xx, yy = np.meshgrid(y, x) # Not typo
        zz = np.arctan2(xx, yy)+np.pi

        bbox0 = Bbox.from_bounds(x0, y0, width, height)
        bbox = TransformedBbox(bbox0, handlebox.get_transform())

        img = BboxImage(bbox, cmap=cmap, norm=norm, data=zz)
        arrow = mpl.patches.FancyArrowPatch((x0, y0+height/2), (x0+width, y0+height/2),
                                 mutation_scale=50, transform=handlebox.get_transform())
        img.set_clip_path(arrow)

        handlebox.add_artist(img)
        return img

plt.legend(loc='best', fontsize=16, framealpha=1, handler_map={PolyCollection: HandlerImage()})

if animate:
    print('Animating plot...')
    ani = animation.FuncAnimation(plt.gcf(), func=update, frames=range(frames), interval=frame_delay)

    if save_image:
        print('Saving animation...')
        ani.save(filename=out_image_fname, writer="pillow", fps=1000/frame_delay)
        # FFwriter = animation.FFMpegWriter(fps=10)
        # ani.save(filename="na-macrostrat.mov", writer=FFwriter)

print('Processing complete!')
if show_plot:
    plt.show()