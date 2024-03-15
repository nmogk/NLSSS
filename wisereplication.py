import requests
import more_itertools
import itertools
from collections import deque
from tqdm import tqdm
import spatialite as sqlite3
import pickle
import csv
from sql_statements import *
from paleobiodb_interface import *

# Settings for 
search_lvl = 5
threshold_distance_deg = 2
taxon_level = 'species' # species, genus, family
env_type = None # None, terr, marine
taxa_filt = None # plantae, prokaryota,eukaryota^plantae
count_global_crossings = True
find_gappers = True # Include taxa which straddle a boundary with any number of series gaps

# Provide the filename for the CSV file
csv_filename = 'output.csv'

column_filename = 'column.pkl'
database_filename = 'paleobiodb.sqlite'

def taxon_field_picker(level):
    if level=='species':
        return rv.SPECIES
    if level=='genus':
        return rv.GENUS
    if level=='family':
        return rv.FAMILY
taxon_field = taxon_field_picker(taxon_level)

# Initialize queries from settings fields
init_sql_statements(taxon_field, threshold_distance_deg)
init_paleobiodb_queries(taxon_level, env_type, taxa_filt)

# Select result labels based on the selected taxon analysis level
total_res_label = 'total_' + taxon_level
label_temp = list('nlsss')
label_temp[2] = list('orpes')[search_lvl-1] # eOnotherm, eRathem, Period (system), Epoch (series), Stage
label_temp[4] = taxon_level[0]

global_temp = label_temp.copy()
global_temp[1] = 'g'

local_label = ''.join(label_temp)
global_label = ''.join(global_temp)

label_temp[-2] = 'j' # j for "Jumping"
local_gap_label = ''.join(label_temp)
global_temp[-2] = 'j'
global_gap_label = ''.join(global_temp)


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

        res = requests.get(api_base+interval_request.format(interval[rv.LEVEL] + 1)+column_parent_fragment.format(interval[rv.MIN_MA], interval[rv.MAX_MA]))
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
    load_db_from_file(conn, database_filename)
    print('spatialite version: ' + conn.execute(spatialite_query).fetchone()[0])

    # Perform spatial queries using SpatiaLite functions
    cursor = conn.cursor()
    
    
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

    save_db_to_file(conn, database_filename)

    id = 1
    result = {}

    print('Processing boundaries...')
    for below, above in tqdm(more_itertools.windowed(column, 2), total=len(column)-1):
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



        cursor.execute(countQuery.format(lowertable))
        res[total_res_label] = cursor.fetchone()[0]

        if find_gappers:

            cursor.executescript(create_union_view(lowertable+'_olderview', [tableName(age[rv.NAME]) for age in itertools.islice(column, 0, id)]))
            conn.commit()
            cursor.executescript(create_union_view(uppertable+'_youngerview', [tableName(age[rv.NAME]) for age in itertools.islice(column, id, None)]))
            conn.commit()

            cursor.execute(copyQuery.format(newtable=lowertable + '_localgappers', table1=lowertable + '_olderview', table2=uppertable + '_youngerview'))
            conn.commit()

            cursor.execute(countQuery.format(lowertable + '_localgappers'))
            res[local_gap_label] = cursor.fetchone()[0]

            if count_global_crossings:
                cursor.execute(copyGlobalQuery.format(newtable=lowertable + '_globalgappers', table1=lowertable + '_olderview', table2=uppertable + '_youngerview'))
                conn.commit()

                cursor.execute(countQuery.format(lowertable + '_globalgappers'))
                res[global_gap_label] = cursor.fetchone()[0]

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
            unionresult = countUnion.format(table1=lowertable, table2=uppertable)
            cursor.execute(countQuery.format(unionresult))
            denom = cursor.fetchone()[0]
            result[id-1][local_label + '_pct'] = 0 if denom == 0 else result[id-1][local_label]/denom
            
            if count_global_crossings:
                result[id-1][global_label + '_pct'] = 0 if denom == 0 else result[id-1][global_label]/denom
            
            if find_gappers:
                unionresult = countUnion.format(table1=lowertable + '_olderview', table2=uppertable + '_youngerview')
                cursor.execute(countQuery.format(unionresult))
                denom = cursor.fetchone()[0]
                result[id-1][local_gap_label + '_pct'] = 0 if denom == 0 else result[id-1][local_gap_label]/denom

                if count_global_crossings:
                    result[id-1][global_gap_label + '_pct'] = 0 if denom == 0 else result[id-1][global_gap_label]/denom

        result[id] = res
        id += 1
    
    save_db_to_file(conn, database_filename)

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

def clearProcessedBoundaries(local=True, glob=True, gappers=True):
    with sqlite3.connect(database_filename) as conn:
        cursor = conn.cursor()
        if local:
            (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_localcrossings')) for interval in column)
        if glob:
            (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_globalcrossings')) for interval in column)
        if gappers:
            (cursor.execute(dropViewQuery.format(tableName(interval[rv.NAME])+'_youngerview')) for interval in column)
            (cursor.execute(dropViewQuery.format(tableName(interval[rv.NAME])+'_olderview')) for interval in column)
            (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_localgappers')) for interval in column)
            if glob:
                (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_globalgappers')) for interval in column)
        conn.commit()

# Export the dictionary of dictionaries to a CSV file
export_dict_of_dicts_to_csv(result, csv_filename)
print(f'Results written to: {csv_filename}')
# for line in column:
#     print(line)