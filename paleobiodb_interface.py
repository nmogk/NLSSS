from strenum import StrEnum
import requests as req

api_base = 'https://paleobiodb.org/data1.2/'

# Translation between meaning and specifid field id returned (or accepted) in the paleobiodb API
rv = StrEnum('ResponseVocab', [('ID', 'oid'), ('NAME', 'nam'), ('MAX_MA', 'eag'), 
                               ('MIN_MA', 'lag'), ('PARENT', 'pid'), ('LEVEL', 'itp'), 
                               ('LAT', 'lat'), ('LON', 'lng'), ('SPECIES', 'tna'), 
                               ('PRECISION', 'prc'), ('FAMILY', 'fml'), ('GENUS', 'gnl'), 
                               ('ENVIRONMENT', 'envtype'), ('FILTER_TAXA', 'base_name'), ('TAXON_ID', 'tid')])

int_vocab = StrEnum('IntervalVocab', [('ID', 'oid'), ('NAME', 'nam'), ('PARENT', 'pid'), ('MAX_MA', 'eag'), ('MIN_MA', 'lag'), ('LEVEL', 'itp')])

clad_vocab = StrEnum('CladogramVocab', [('ID', 'oid'), ('NAME', 'nam'), ('PARENT_ID', 'par'), ('FLAG', 'flg')])

interval_request = 'intervals/list.json?scale=1'
column_parent_fragment = '&min_ma={}&max_ma={}'

occurrence_request = ''

def init_paleobiodb_queries(taxon_level, env_type=None, taxa_filt=None):
    global occurrence_request
    occurrence_request = ('occs/list.json?interval_id={}&pres=regular&show=acconly,class,coords,loc&idreso=' + taxon_level + 
            ('&' + rv.ENVIRONMENT + '=' + env_type if env_type is not None else '') + 
            ('&' + rv.FILTER_TAXA + '=' + taxa_filt if taxa_filt is not None else ''))


def query_occs_by_taxon(taxa_filt):
    taxon_occs_request = f'occs/list.json?{rv.FILTER_TAXA}={taxa_filt}&show=acconly,coords'

    res = req.get(api_base + taxon_occs_request)
    data = res.json()['records']
    return data

def query_geological_intervals(min_ma=None, max_ma=None):
    query_extra = column_parent_fragment.format(min_ma, max_ma) if min_ma is not None and max_ma is not None else ''

    res = req.get(api_base + interval_request + query_extra)
    data = res.json()['records']
    return data

def query_cladogram(taxon_name):
    taxonomy_request = 'occs/taxa.json?base_name={}'

    res = req.get(api_base + taxonomy_request.format(taxon_name))
    data = res.json()['records']
    return data