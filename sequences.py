import pickle
from paleobiodb_interface import rv
import requests as req
import json

baseUrl = 'https://macrostrat.org/api/'

column_filename = 'column.pkl'

with open(column_filename, 'rb') as f:
    column = pickle.load(f)

    for interval in column:
        res = req.get(baseUrl+f'sections?interval_name={interval[rv.NAME]}')
        try:
            packages = len(res.json()['success']['data'])
            print(f'{interval[rv.NAME]}, {packages}')
        except:
            print(f'{interval[rv.NAME]}, -')