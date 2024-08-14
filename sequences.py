import pickle
from paleobiodb_interface import rv
import requests as req
import json

baseUrl = 'https://macrostrat.org/api/'

column_filename = 'column.pkl'

with open(column_filename, 'rb') as f:
    column = pickle.load(f)

    for interval in range(0, 540, 6):
        # res = req.get(baseUrl+f'sections?interval_name={interval[rv.NAME]}&environ_class=marine')
        res = req.get(baseUrl+f'sections?age_bottom={interval}&age_top={interval+10}')
        try:
            # Filter out 0 thickness packages
            filtered = [x for x in res.json()['success']['data'] if x['max_thick'] != '0.00']
            packages = len(filtered)
            print(f'{interval}, {packages}')
        except:
            print(f'{interval}, -')