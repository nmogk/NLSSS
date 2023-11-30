import requests
import more_itertools
from collections import deque
from strenum import StrEnum
from tqdm import tqdm
# import spatialite as sqlite3
import sqlite3
import pickle

search_lvl = 5
threshold_distance_km = 100
taxon_level = 'species' # species, genus, family
env_type = None # None, terr, marine envtype

rv = StrEnum('ResponseVocab', [('ID', 'oid'), ('NAME', 'nam'), ('MAX_MA', 'eag'), ('MIN_MA', 'lag'), ('PARENT', 'pid'), ('LEVEL', 'lvl'), ('LAT', 'lat'), ('LON', 'lng'), ('SPECIES', 'tna'), ('PRECISION', 'prc'), ('FAMILY', 'fml'), ('GENUS', 'gnl')])

column_filename = 'column.pkl'
database_filename = 'paleobiodb.sqlite'
api_base = 'https://paleobiodb.org/data1.2/'

# Download files for all species in the lower and upper intervals of a boundary
# Delete all records identified less precicely than species level
# Consider at the species level those identified to subspecies level
# Delete all records located less precicely than interval
# Accept local stratigraphic units where boundaries are close to global boundaries (2My for most, 1My for units 2 My and shorter, 0My for top 2)
# Combine files
# Delete species found in only one side of boundary
# Delete species (sic, occurrences?) located more than 2 degrees away from any other specimen on the other side of the boundary
# Count species

def queryColumn():
    interval_request = 'intervals/list.json?scale=1&level={}'
    parent_fragment = '&min_ma={}&max_ma={}'

    # Initial query to get highest level intervals
    res = requests.get(api_base+interval_request.format(1))
    seedData = res.json()

    # Load intervals into stack LIFO (oldest on top)
    stack = deque()
    for record in seedData['records']:
        stack.append(record)

    def checkSubintervals(parent, childList):
        '''childList must be sorted youngest to oldest'''
        pointer = parent[rv.MIN_MA]
        for child in childList:
            if child[rv.PARENT] != parent[rv.ID]:
                continue
            if child[rv.MIN_MA] != pointer:
                return False
            pointer = child[rv.MAX_MA]
        return parent[rv.MAX_MA] == pointer

    column = deque()
    t = tqdm(total=115)
    while True:
        if len(stack) <= 0:
            break
        interval = stack.pop()

        if interval[rv.LEVEL] >= search_lvl:
            column.append(interval)
            t.update(1)
            continue

        res = requests.get(api_base+interval_request.format(interval[rv.LEVEL] + 1)+parent_fragment.format(interval[rv.MIN_MA], interval[rv.MAX_MA]))
        subintervals = res.json()

        if checkSubintervals(interval, subintervals['records']):
            for subint in subintervals['records']:
                stack.append(subint)
        else:
            column.append(interval)
            t.update(1)
    t.close()
    return column

print('Loading column information...')
try:
    with open(column_filename, 'rb') as f:
        column = pickle.load(f)
except FileNotFoundError as err:
    column = queryColumn()
    with open(column_filename, 'wb') as f:
        pickle.dump(column, f)

# Connect to a SQLite database (which includes SpatiaLite)
with sqlite3.connect(database_filename) as conn:
    # print(conn.execute('SELECT spatialite_version()').fetchone()[0])

    # Enable the SpatiaLite extension
    conn.enable_load_extension(True)
    # conn.load_extension('mod_spatialite')
    
    # Perform spatial queries using SpatiaLite functions
    cursor = conn.cursor()
    
    occurrence_request = 'occs/list.json?interval_id={}&pres=regular&show=acconly,class,coords,loc&idreso=' + taxon_level
    check_table_query = 'SELECT 1 FROM sqlite_master WHERE type="table" AND name="{}"'
    create_table_query = 'CREATE TABLE {}(' + ', '.join((rv.ID, rv.LAT, rv.LON, rv.PRECISION, rv.SPECIES, rv.GENUS, rv.FAMILY)) +  ')'
    insert_query = 'INSERT INTO {} VALUES (?, ?, ?, ?, ?, ?, ?)'
    
    def get_insert_values(occurrence):
        return (occurrence[rv.ID] , 
                  float(occurrence[rv.LAT]), 
                  float(occurrence[rv.LON]), 
                  occurrence[rv.PRECISION], 
                  occurrence[rv.SPECIES], 
                  occurrence[rv.GENUS], 
                  occurrence[rv.FAMILY])
    
    print('Downloading fossil occurrence data...')
    for interval in tqdm(column):
        tablename = interval[rv.NAME].replace(' ', '_')

        cursor.execute(check_table_query.format(tablename))
        if cursor.fetchone() is not None:
            continue
    
        res = requests.get(api_base+occurrence_request.format(interval[rv.ID]))
        occs = res.json()['records']
        # Load result into database
        cursor.execute(create_table_query.format(tablename))
        cursor.executemany(insert_query.format(tablename), (get_insert_values(occ) for occ in occs))
        conn.commit()

print('Processing boundaries...')
for below, above in tqdm(more_itertools.windowed(column, 2), total=114):
    pass
    # Algorithm:
    # Count total unique species
    # Count unique species which globally cross boundary
    # Starting from bottom, proceding upwards:
    #     Delete occurrences of species unique to lower unit
    #     Delete occurrences in lower unit farther than threshold distance from occurrences of the same species in upper unit
    #     Count remaining unique species
    # Save results data
    # {id: {bdry_name:, total_species:, ngsss:, nlsss:}}
    # Export to csv

# query = f"""
#     SELECT *
#     FROM your_table
#     WHERE ST_Distance(
#         MakePoint({point_longitude}, {point_latitude}, 4326),
#         geom_column
#     ) <= {threshold_distance_km} * 1000
# """

# cursor.execute(query)

# result = cursor.fetchall()

# for row in result:
#     print(row)  # Handle nearby points

# Close the connection
# conn.close()

# for line in column:
#     print(line)