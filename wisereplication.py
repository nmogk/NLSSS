import requests
import more_itertools
import itertools
from collections import deque
from strenum import StrEnum
from tqdm import tqdm
import spatialite as sqlite3
import pickle
import csv

search_lvl = 5
threshold_distance_deg = 2
taxon_level = 'species' # species, genus, family
env_type = None # None, terr, marine
taxa_filt = None # plantae, prokaryota,eukaryota^plantae
count_global_crossings = True
find_gappers = True # Include taxa which straddle a boundary with any number of series gaps

# Provide the filename for the CSV file
csv_filename = 'output.csv'

rv = StrEnum('ResponseVocab', [('ID', 'oid'), ('NAME', 'nam'), ('MAX_MA', 'eag'), ('MIN_MA', 'lag'), ('PARENT', 'pid'), ('LEVEL', 'lvl'), ('LAT', 'lat'), ('LON', 'lng'), ('SPECIES', 'tna'), ('PRECISION', 'prc'), ('FAMILY', 'fml'), ('GENUS', 'gnl'), ('ENVIRONMENT', 'envtype'), ('FILTER_TAXA', 'base_name')])

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

def create_union_view(view_name, table_names):
    query = f"DROP VIEW IF EXISTS {view_name};\n"  # Drop view if it exists
    query += f"CREATE VIEW {view_name} AS\n"
    for i, table_name in enumerate(table_names):
        if i > 0:
            query += "UNION ALL\n"
        query += f"SELECT * FROM {table_name}\n"
    return query

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
    
    occurrence_request = ('occs/list.json?interval_id={}&pres=regular&show=acconly,class,coords,loc&idreso=' + taxon_level + 
            ('&' + rv.ENVIRONMENT + '=' + env_type if env_type is not None else '') + 
            ('&' + rv.FILTER_TAXA + '=' + taxa_filt if taxa_filt is not None else ''))
    check_table_query = 'SELECT 1 FROM sqlite_schema WHERE type="table" AND name="{}"'
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

    # Save in-memory database to disk
    # Connect to the SQLite database on disk
    conn_disk = sqlite3.connect(database_filename)
    # Backup the data from memory to disk
    conn.backup(conn_disk)
    # Close the database connections
    conn_disk.close()

    copyQuery = (
    'CREATE TABLE IF NOT EXISTS {newtable} AS ' +
    'SELECT * FROM {table1} ' +
    'WHERE EXISTS (' +
        'SELECT 1 ' +
        'FROM {table2} ' +
        'WHERE {table1}.' + taxon_field + ' = {table2}.' + taxon_field +
        ' AND ST_Distance({table1}.location, {table2}.location) <= ' + str(threshold_distance_deg) + ' )'
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
    countUnion = '(SELECT ' + taxon_field + ' FROM {table1} UNION SELECT ' + taxon_field + ' FROM {table2})'

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

        # Select result labels based on the selected taxon analysis level
        total_res_label = 'total_' + taxon_level
        label_temp = list('nlsss')
        label_temp[4] = taxon_level[0]
        local_label = ''.join(label_temp)
        label_temp[1] = 'g'
        global_label = ''.join(label_temp)

        cursor.execute(countQuery.format(lowertable))
        res[total_res_label] = cursor.fetchone()[0]

        if find_gappers:
            label_temp = list(local_label)
            label_temp[-2] = 'j' # j for "Jumping"
            local_gap_label = ''.join(label_temp)

            cursor.executescript(create_union_view(lowertable+'_olderview', [tableName(age[rv.NAME]) for age in itertools.islice(column, 0, id)]))
            conn.commit()
            cursor.executescript(create_union_view(uppertable+'_youngerview', [tableName(age[rv.NAME]) for age in itertools.islice(column, id, None)]))
            conn.commit()

            cursor.execute(copyQuery.format(newtable=lowertable + '_localgappers', table1=lowertable + '_olderview', table2=uppertable + '_youngerview'))
            conn.commit()

            cursor.execute(countQuery.format(lowertable + '_localgappers'))
            res[local_gap_label] = cursor.fetchone()[0]

            if count_global_crossings:
                label_temp[1] = 'g'
                global_gap_label = ''.join(label_temp)
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

dropTableQuery = 'DROP TABLE IF EXISTS {}'
def clearProcessedBoundaries(local=True, glob=True):
    with sqlite3.connect(database_filename) as conn:
        cursor = conn.cursor()
        if local:
            (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_localcrossings')) for interval in column)
        if glob:
            (cursor.execute(dropTableQuery.format(tableName(interval[rv.NAME])+'_globalcrossings')) for interval in column)
        conn.commit()

# Export the dictionary of dictionaries to a CSV file
export_dict_of_dicts_to_csv(result, csv_filename)
print(f'Results written to: {csv_filename}')
# for line in column:
#     print(line)