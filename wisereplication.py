import requests
import more_itertools
import itertools
from collections import deque
from tqdm import tqdm
import spatialite as sqlite3
import pickle
import csv
import sql_statements as sql
import paleobiodb_interface as pbdb
from paleobiodb_interface import rv
from multiprocess import Manager
from strenum import StrEnum
from enum import auto
import sys

class TimeLevel(StrEnum):
    eon = auto()
    era = auto()
    period = auto()
    epoch = auto()
    # subepoch = 5
    age = auto()
    # subage = 7
    # zone = 8

    def index(self):
        cls = self.__class__
        members = list(cls)
        return members.index(self) + 1

    def next(self):
        cls = self.__class__
        members = list(cls)
        index = members.index(self) + 1
        if index >= len(members):
            raise StopIteration('end of enumeration reached')
        return members[index]

    @classmethod
    def abbreviate_levels(cls, level):
        if level == cls.eon:
            return 'o'
        if level == cls.era:
            return 'r'
        if level == cls.period:
            return 'p'
        if level == cls.epoch:
            return 'e'
        if level == cls.age:
            return 's'

# Settings for 
search_lvl = TimeLevel.age # How many levels deep to generate column names: 1-eon, 5-stage
threshold_distance_deg = 2
taxon_level = 'species' # species, genus, family
env_type = None # None, terr, marine
taxa_filt = None # plantae, prokaryota,eukaryota^plantae
count_global_crossings = True
find_gappers = False # Include taxa which straddle a boundary with any number of series gaps

# Provide the filename for the CSV file
csv_filename = 'fbwg_nlsss_base.csv'

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
sql.init_sql_statements(taxon_field, threshold_distance_deg)
pbdb.init_paleobiodb_queries(taxon_level, env_type, taxa_filt)

# Select result labels based on the selected taxon analysis level
total_res_label = 'total_' + taxon_level
label_temp = list('nlsss')
label_temp[2] = TimeLevel.abbreviate_levels(search_lvl)
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
    res = requests.get(pbdb.api_base+pbdb.interval_request)
    seedData = res.json()

    # Load intervals into stack LIFO (oldest on top)
    stack = deque()
    for record in seedData['records']:
        if TimeLevel[record[rv.LEVEL]] == TimeLevel.eon:
            stack.append(record)

    def checkSubintervals(parent, childList):
        '''childList must be sorted youngest to oldest'''
        pointer = parent[rv.MIN_MA]
        for child in childList:
            if TimeLevel[child[rv.LEVEL]].index() != TimeLevel(parent[rv.LEVEL]).next().index() or child[rv.PARENT] != parent[rv.ID]:
                continue
            if child[rv.MIN_MA] != pointer:
                return False
            pointer = child[rv.MAX_MA]
        return parent[rv.MAX_MA] == pointer

    column = deque()
    t = tqdm(total=117)
    while True:
        if len(stack) <= 0:
            break
        interval = stack.pop()

        if TimeLevel[interval[rv.LEVEL]].index() >= search_lvl.index():
            column.append(interval)
            t.update(1)
            continue

        res = requests.get(pbdb.api_base+pbdb.interval_request+pbdb.column_parent_fragment.format(interval[rv.MIN_MA], interval[rv.MAX_MA]))
        subintervals = res.json()

        if checkSubintervals(interval, subintervals['records']):
            for subint in subintervals['records']:
                if TimeLevel[subint[rv.LEVEL]].index() == TimeLevel(interval[rv.LEVEL]).next().index():
                    stack.append(subint)
        else:
            column.append(interval)
            t.update(1)
    t.close()
    return column

def tableName(textname):
    return textname.replace(' ', '_').lower()

def retreive_paleobiodb_data(column):
    # Connect to a SQLite database (which includes SpatiaLite)
    with sqlite3.connect(':memory:') as conn:
        success = True

        # Attach the database from disk to memory. For a variety of setups, this will shave minutes off processing time
        sql.load_db_from_file(conn, database_filename)
        print('spatialite version: ' + conn.execute(sql.spatialite_query).fetchone()[0])

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

            cursor.execute(sql.check_table_query.format(tablename))
            if cursor.fetchone() is not None:
                continue
            
            res = requests.get(pbdb.api_base+pbdb.occurrence_request.format(interval[rv.ID]))
            try:
                occs = res.json()['records']
            except KeyError:
                print(f'Error returned when querying PaleoBioDB for {interval[rv.NAME]}. Please refresh geological column and download data again.')
                success = False
                continue
            # Load result into database
            cursor.execute(sql.create_table_query.format(tablename))
            cursor.executemany(sql.insert_query.format(tablename), (get_insert_values(occ) for occ in occs))
            conn.commit()

        sql.save_db_to_file(conn, database_filename)
        return success

def find_bounary_crossers(column):
    print('Processing boundaries...')
    manager = Manager()
    result = manager.dict()
    
    def worker_tasks(input):
        # Algorithm:
        # Count total unique species
        # Count unique species which globally cross boundary
        # Starting from bottom, proceding upwards:
        #     Create new table of occurrences of species which cross the boundary
        #     Create new table of occurrences in lower unit which cross the boundary and are closer than threshold distance from occurrences of the same species in upper unit
        #     Count remaining unique species globally and locally
        # Save results data
        # {id: {boundary: (name), total_species:, ngsss:, nlsss:, ngsjs:, nlsjs:, ngsss_pct:, nlsss_pct:, ngsjs_pct:, nlsjs_pct:}}
        id, window = input
        below, above = window
        with sqlite3.connect(database_filename, 60) as conn:
            conn.execute("PRAGMA journal_mode=WAL")
            # Perform spatial queries using SpatiaLite functions
            cursor = conn.cursor()
            cursor.execute("BEGIN IMMEDIATE")

            lowertable = tableName(below[rv.NAME])
            uppertable = tableName(above[rv.NAME])
            res = {'boundary': '/'.join((below[rv.NAME], above[rv.NAME]))}

            cursor.execute(sql.countQuery.format(lowertable))
            res[total_res_label] = cursor.fetchone()[0]

            if find_gappers:
                cursor.executescript(sql.create_union_view(lowertable+'_olderview', [tableName(age[rv.NAME]) for age in itertools.islice(column, 0, id)]))
                conn.commit()
                cursor.executescript(sql.create_union_view(uppertable+'_youngerview', [tableName(age[rv.NAME]) for age in itertools.islice(column, id, None)]))
                conn.commit()

                cursor.execute(sql.copyQuery.format(newtable=lowertable + '_localgappers', table1=lowertable + '_olderview', table2=uppertable + '_youngerview'))
                conn.commit()

                cursor.execute(sql.countQuery.format(lowertable + '_localgappers'))
                res[local_gap_label] = cursor.fetchone()[0]

                if count_global_crossings:
                    cursor.execute(sql.copyGlobalQuery.format(newtable=lowertable + '_globalgappers', table1=lowertable + '_olderview', table2=uppertable + '_youngerview'))
                    conn.commit()

                    cursor.execute(sql.countQuery.format(lowertable + '_globalgappers'))
                    res[global_gap_label] = cursor.fetchone()[0]

            if count_global_crossings:
                # This query only deletes occurrences without any members that cross the boundary
                # Only needed if we want to count global crossings, since the distance query will also delete occurrences without any crossings
                cursor.execute(sql.copyGlobalQuery.format(newtable=lowertable + '_globalcrossings', table1=lowertable, table2=uppertable))
                conn.commit()

                cursor.execute(sql.countQuery.format(lowertable + '_globalcrossings'))
                res[global_label] = cursor.fetchone()[0]

            # Delete occurrences of species unique to lower unit or without members above the boundary closer than the threshold distance.
            cursor.execute(sql.copyQuery.format(newtable=lowertable + '_localcrossings', table1=lowertable, table2=uppertable))
            conn.commit()

            cursor.execute(sql.countQuery.format(lowertable + '_localcrossings'))
            res[local_label] = cursor.fetchone()[0]

            result[id] = res

    ppool = manager.Pool(1)
    with tqdm(total=len(column)-1) as pbar:
        for _ in ppool.imap_unordered(worker_tasks, list(enumerate(more_itertools.windowed(column, 2), 1))):
            pbar.update()

    ppool.close()
    ppool.join()

    # with sqlite3.connect(':memory:') as conn:
    #     sql.save_db_to_file(conn, database_filename)

    return result

def overlap_statistics(column, result):
    with sqlite3.connect(database_filename) as conn:
        # Perform spatial queries using SpatiaLite functions
        cursor = conn.cursor()

        for id, (below, above) in tqdm(enumerate(more_itertools.windowed(column, 2)), total=len(column)-1):
            lowertable = tableName(below[rv.NAME])
            uppertable = tableName(above[rv.NAME])
            if id > 0:
                unionresult = sql.countUnion.format(table1=lowertable, table2=uppertable)
                cursor.execute(sql.countQuery.format(unionresult))
                denom = cursor.fetchone()[0]
                result[id][local_label + '_pct'] = 0 if denom == 0 else result[id][local_label]/denom

                if count_global_crossings:
                    result[id][global_label + '_pct'] = 0 if denom == 0 else result[id][global_label]/denom

                if find_gappers:
                    unionresult = sql.countUnion.format(table1=lowertable + '_olderview', table2=uppertable + '_youngerview')
                    cursor.execute(sql.countQuery.format(unionresult))
                    denom = cursor.fetchone()[0]
                    result[id][local_gap_label + '_pct'] = 0 if denom == 0 else result[id][local_gap_label]/denom

                    if count_global_crossings:
                        result[id][global_gap_label + '_pct'] = 0 if denom == 0 else result[id][global_gap_label]/denom

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
            (cursor.execute(sql.dropTableQuery.format(tableName(interval[rv.NAME])+'_localcrossings')) for interval in column)
        if glob:
            (cursor.execute(sql.dropTableQuery.format(tableName(interval[rv.NAME])+'_globalcrossings')) for interval in column)
        if gappers:
            (cursor.execute(sql.dropViewQuery.format(tableName(interval[rv.NAME])+'_youngerview')) for interval in column)
            (cursor.execute(sql.dropViewQuery.format(tableName(interval[rv.NAME])+'_olderview')) for interval in column)
            (cursor.execute(sql.dropTableQuery.format(tableName(interval[rv.NAME])+'_localgappers')) for interval in column)
            if glob:
                (cursor.execute(sql.dropTableQuery.format(tableName(interval[rv.NAME])+'_globalgappers')) for interval in column)
        conn.commit()

if __name__ == "__main__":
    print('Loading column information...')
    try:
        with open(column_filename, 'rb') as f:
            column = pickle.load(f)
    except FileNotFoundError as err:
        column = queryColumn()
        with open(column_filename, 'wb') as f:
            pickle.dump(column, f)
    # for line in column:
    #     print(line)

    download_success = retreive_paleobiodb_data(column) # Single process
    if not download_success:
        sys.exit()

    result = find_bounary_crossers(column) # Multiprocess
    overlap_statistics(column, result) # Single process

    # Export the dictionary of dictionaries to a CSV file
    export_dict_of_dicts_to_csv(result, csv_filename)
    print(f'Results written to: {csv_filename}')