import requests
import more_itertools
from collections import deque
from strenum import StrEnum
from tqdm import tqdm
import spatialite as sqlite3
import pickle
import csv

search_lvl = 5
threshold_distance_km = 220
taxon_level = 'species' # species, genus, family
env_type = None # None, terr, marine
count_global_crossings = True

# Provide the filename for the CSV file
csv_filename = 'output.csv'

rv = StrEnum('ResponseVocab', [('ID', 'oid'), ('NAME', 'nam'), ('MAX_MA', 'eag'), ('MIN_MA', 'lag'), ('PARENT', 'pid'), ('LEVEL', 'lvl'), ('LAT', 'lat'), ('LON', 'lng'), ('SPECIES', 'tna'), ('PRECISION', 'prc'), ('FAMILY', 'fml'), ('GENUS', 'gnl'), ('ENVIRONMENT', 'envtype')])

column_filename = 'column.pkl'
database_filename = 'paleobiodb.sqlite'
api_base = 'https://paleobiodb.org/data1.2/'

def taxon_field_picker(level):
    match level:
        case 'species':
            return rv.SPECIES
        case 'genus':
            return rv.GENUS
        case 'family':
            return rv.FAMILY
taxon_field = taxon_field_picker(taxon_level)

# Original Wise algorithm
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

def tableName(textname):
    return textname.replace(' ', '_').lower()

# Connect to a SQLite database (which includes SpatiaLite)
with sqlite3.connect(':memory:') as conn:

    # Attach the database from disk to memory. For a variety of setups, this will shave minutes off processing time
    # Connect to the SQLite database on disk
    conn_disk = sqlite3.connect(database_filename)
    # Backup the data from disk to memory
    conn_disk.backup(conn)
    # Close the database connections
    conn_disk.close()
    print('spatialite version: ' + conn.execute('SELECT spatialite_version()').fetchone()[0])

    # Perform spatial queries using SpatiaLite functions
    cursor = conn.cursor()
    
    occurrence_request = 'occs/list.json?interval_id={}&pres=regular&show=acconly,class,coords,loc&idreso=' + taxon_level + ('&' + rv.ENVIRONMENT + '=' + env_type if env_type is not None else '')
    check_table_query = 'SELECT 1 FROM sqlite_master WHERE type="table" AND name="{}"'
    create_table_query = 'CREATE TABLE {}(' + ', '.join((rv.ID, 'location', rv.PRECISION, rv.SPECIES, rv.GENUS, rv.FAMILY)) +  ')'
    insert_query = 'INSERT INTO {} VALUES (?, MakePoint(? ,? ,4326), ?, ?, ?, ?)'
    
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
        tablename = tableName(interval[rv.NAME])

        cursor.execute(check_table_query.format(tablename))
        if cursor.fetchone() is not None:
            continue
    
        res = requests.get(api_base+occurrence_request.format(interval[rv.ID]))
        occs = res.json()['records']
        # Load result into database
        cursor.execute(create_table_query.format(tablename))
        cursor.executemany(insert_query.format(tablename), (get_insert_values(occ) for occ in occs))
        conn.commit()

    copyQuery = (
    'CREATE TABLE IF NOT EXISTS {newtable} AS ' +
    'SELECT * FROM {table1} ' +
    'WHERE EXISTS (' +
        'SELECT 1 ' +
        'FROM {table2} ' +
        'WHERE {table1}.' + taxon_field + ' = {table2}.' + taxon_field +
        ' AND ST_Distance({table1}.location, {table2}.location) <= ' + str(threshold_distance_km) + ' * 1000)'
    )

    copyGlobalQuery = (
    'CREATE TABLE IF NOT EXISTS {newtable} AS ' +
    'SELECT * FROM {table1} ' +
    'WHERE EXISTS (' +
        'SELECT 1 ' +
        'FROM {table2} ' +
        'WHERE {table1}.' + taxon_field + ' = {table2}.' + taxon_field + ')'
    )

    countQuery = 'SELECT COUNT(DISTINCT ' + taxon_field + ') FROM {}'

    id = 1
    result = {}

    print('Processing boundaries...')
    for below, above in tqdm(more_itertools.windowed(column, 2), total=114):
        # Algorithm:
        # Count total unique species
        # Count unique species which globally cross boundary
        # Starting from bottom, proceding upwards:
        #     Create new table of occurrences of species which cross the boundary
        #     Create new table of occurrences in lower unit which cross the boundary and are closer than threshold distance from occurrences of the same species in upper unit
        #     Count remaining unique species globally and locally
        # Save results data
        # {id: {bdry_name:, total_species:, ngsss:, nlsss:}}
        # Export to csv
        lowertable = tableName(below[rv.NAME])
        uppertable = tableName(above[rv.NAME])
        res = {'boundary': '/'.join((below[rv.NAME], above[rv.NAME]))}

        # Select result labels based on the selected taxon analysis level
        total_res_label = 'total_' + taxon_level
        label_temp = list('nlsss')
        label_temp[4] = taxon_level[0]
        local_label = ''.join(label_temp)
        label_temp[1] = 'g'
        global_label = ''.join(label_temp)

        cursor.execute(countQuery.format(lowertable))
        res[total_res_label] = cursor.fetchone()[0]

        if count_global_crossings:
            # This query only deletes occurrences without any members that cross the boundary
            # Only needed if we want to count global crossings, since the distance query will also delete occurrences without any crossings
            cursor.execute(copyGlobalQuery.format(newtable=lowertable + '_globalcrossings', table1=lowertable, table2=uppertable))
            conn.commit()

            cursor.execute(countQuery.format(lowertable + '_globalcrossings'))
            res[global_label] = cursor.fetchone()[0]

        # Delete occurrences of species unique to lower unit or without members above the boundary closer than the threshold distance.
        cursor.execute(copyQuery.format(newtable=lowertable + '_localcrossings', table1=lowertable, table2=uppertable))
        conn.commit()

        cursor.execute(countQuery.format(lowertable + '_localcrossings'))
        res[local_label] = cursor.fetchone()[0]

        if id > 1:
            denom = result[id-1][total_res_label] + res[total_res_label]
            result[id-1][local_label + '_pct'] = 0 if denom == 0 else result[id-1][local_label]/denom
            if count_global_crossings:
                result[id-1][global_label + '_pct'] = 0 if denom == 0 else result[id-1][global_label]/denom

        result[id] = res
        id += 1
    
    # Detach the in-memory database and save its changes to disk
    # Connect to the SQLite database on disk
    conn_disk = sqlite3.connect(database_filename)
    # Backup the data from memory to disk
    conn.backup(conn_disk)
    # Close the database connections
    conn_disk.close()

def export_dict_of_dicts_to_csv(data, csv_filename):
    # Extract headers from the first dictionary
    headers = list(data[next(iter(data))].keys())

    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['bdry_no'] + headers)
        
        writer.writeheader()

        for key, inner_dict in data.items():
            row = {'bdry_no': key}
            row.update(inner_dict)
            writer.writerow(row)

# Export the dictionary of dictionaries to a CSV file
export_dict_of_dicts_to_csv(result, csv_filename)
print(f'Results written to: {csv_filename}')
# for line in column:
#     print(line)