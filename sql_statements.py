from paleobiodb_interface import rv
import sqlite3

# Check to see if spatialite is installed in the database
spatialite_query = 'SELECT spatialite_version()' 

# Check to see if a table already exists
check_table_query = 'SELECT 1 FROM sqlite_schema WHERE type="table" AND name="{}"' 

# Table creation query with application-specific field IDs
create_table_query = 'CREATE TABLE {}(' + ', '.join((rv.ID, 'location', rv.PRECISION, rv.SPECIES, rv.GENUS, rv.FAMILY)) +  ')'

# Generic insert consistent with column definition above
insert_query = 'INSERT INTO {} VALUES (?, MakePoint(? ,? ,4326), ?, ?, ?, ?)' 

def create_union_view(view_name, table_names):
    '''Create a view which includes all entries from a list of tables. Drops existing view before creating this one.'''
    query = f"DROP VIEW IF EXISTS {view_name};\n"  # Drop view if it exists
    query += f"CREATE VIEW {view_name} AS\n"
    for i, table_name in enumerate(table_names):
        if i > 0:
            query += "UNION ALL\n"
        query += f"SELECT * FROM {table_name}\n"
    return query

# These queries need initialization
copyQuery, copyGlobalQuery, countQuery, countUnion

def init_sql_statements(taxon_field, threshold_distance_deg):
    '''Initialize statements which require static setting information (specifically, taxon level and spatial search distance) as part of the query'''

    # Create a table that finds occurrences in two specified tables which match in taxon and are within a specified distance of each other
    global copyQuery
    copyQuery = (
    'CREATE TABLE IF NOT EXISTS {newtable} AS ' +
    'SELECT * FROM {table1} ' +
    'WHERE EXISTS (' +
        'SELECT 1 ' +
        'FROM {table2} ' +
        'WHERE {table1}.' + taxon_field + ' = {table2}.' + taxon_field +
        ' AND ST_Distance({table1}.location, {table2}.location) <= ' + str(threshold_distance_deg) + ' )'
    )

    # Create a table that finds occurrences in two specified tables which match in taxon
    global copyGlobalQuery
    copyGlobalQuery = (
    'CREATE TABLE IF NOT EXISTS {newtable} AS ' +
    'SELECT * FROM {table1} ' +
    'WHERE EXISTS (' +
        'SELECT 1 ' +
        'FROM {table2} ' +
        'WHERE {table1}.' + taxon_field + ' = {table2}.' + taxon_field + ')'
    )

    # Count distinct taxa (as opposed to occurrences) in this table
    global countQuery
    countQuery = 'SELECT COUNT(DISTINCT ' + taxon_field + ') FROM {}'

    # Count distinct taxa (as opposed to occurrences) in two tables simultaneously
    global countUnion
    countUnion = '(SELECT ' + taxon_field + ' FROM {table1} UNION SELECT ' + taxon_field + ' FROM {table2})'

# Generic drop table
dropTableQuery = 'DROP TABLE IF EXISTS {}'
dropViewQuery = 'DROP VIEW IF EXISTS {}'

def load_db_from_file(conn, fname):
    # Connect to the SQLite database on disk
    with sqlite3.connect(fname) as conn_disk:
        # Backup the data from disk to memory
        conn_disk.backup(conn)

def save_db_to_file(conn, fname):
    # Save in-memory database to disk
    # Connect to the SQLite database on disk
    with sqlite3.connect(fname) as conn_disk:
        # Backup the data from memory to disk
        conn.backup(conn_disk)
