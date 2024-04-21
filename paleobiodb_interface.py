from strenum import StrEnum

api_base = 'https://paleobiodb.org/data1.2/'

# Translation between meaning and specifid field id returned (or accepted) in the paleobiodb API
rv = StrEnum('ResponseVocab', [('ID', 'oid'), ('NAME', 'nam'), ('MAX_MA', 'eag'), 
                               ('MIN_MA', 'lag'), ('PARENT', 'pid'), ('LEVEL', 'itp'), 
                               ('LAT', 'lat'), ('LON', 'lng'), ('SPECIES', 'tna'), 
                               ('PRECISION', 'prc'), ('FAMILY', 'fml'), ('GENUS', 'gnl'), 
                               ('ENVIRONMENT', 'envtype'), ('FILTER_TAXA', 'base_name')])


interval_request = 'intervals/list.json?scale=1'
column_parent_fragment = '&min_ma={}&max_ma={}'

occurrence_request = ''

def init_paleobiodb_queries(taxon_level, env_type=None, taxa_filt=None):
    global occurrence_request
    occurrence_request = ('occs/list.json?interval_id={}&pres=regular&show=acconly,class,coords,loc&idreso=' + taxon_level + 
            ('&' + rv.ENVIRONMENT + '=' + env_type if env_type is not None else '') + 
            ('&' + rv.FILTER_TAXA + '=' + taxa_filt if taxa_filt is not None else ''))
