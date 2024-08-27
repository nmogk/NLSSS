import pickle
from paleobiodb_interface import rv
import requests as req
import json
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

baseUrl = 'https://macrostrat.org/api/'

data_filename = 'packages_all_Ages.pkl'
use_stages = True
do_smooth = True
step = 6

if use_stages:
    sloss_ms = np.array([467, 419, 324, 132, 60])
    peters_ms = np.array([252, 199])
    mogk_ms = np.array([441, 294, 86, 20])

else:
    sloss_ms = np.array([462, 396, 324, 132, 60])
    peters_ms = np.array([252, 186])
    mogk_ms = np.array([366, 24])

try:
    with open(data_filename, 'rb') as f:
        x, y = pickle.load(f)

except:

    if use_stages:
        res = req.get(f'{baseUrl}defs/intervals?timescale_id=1')
        data = res.json()['success']['data']
        queries = [x['name'] for x in data]
        x = [x['t_age'] for x in data]
    else:
        x = np.arange(0, 540, step)

    y = np.zeros_like(x)
    with open(data_filename, 'wb') as f:
        for i, interval in enumerate(tqdm(x)):
            # res = req.get(baseUrl+f'sections?interval_name={interval[rv.NAME]}&environ_class=marine')
            
            if use_stages:
                res = req.get(f'{baseUrl}sections?interval_name={queries[i]}')#&environ_class=non-marine')
            else:
                res = req.get(baseUrl+f'sections?age_top={interval}&age_bottom={interval+step}')#&environ_class=non-marine')

            # Filter out 0 thickness packages
            filtered = [x for x in res.json()['success']['data'] if x['max_thick'] != '0.00']
            packages = len(filtered)#np.sum([x['col_area'] for x in filtered]) #len(filtered)
            y[i] = packages
        pickle.dump((x,y), f)

if do_smooth:
    window = 5
    kernel = np.array([0.1, .2, 0.4,  .2 ,0.1])#np.ones(window)/window
    y = np.convolve(y, kernel, mode='same')
# print(x)
# print(y)

plt.title(f'Number of marine gap bound packages{", smoothed" if do_smooth  else ""}')
axes = plt.gca()
axes.set_xlim([0, 540])
axes.set_xlabel('Ma')
axes.set_ylabel('Package count')
for bdry in sloss_ms:
    plt.axvline(bdry, color='gray', linestyle='--', label='Sloss (1963) boundaries')
for bdry in peters_ms:
    plt.axvline(bdry, color='blue', linestyle='--', label='Peters (2008) boundaries')
for bdry in mogk_ms:
    plt.axvline(bdry, color='gold', linestyle=':', label='Additional boundaries')
plt.plot(x,y, label=f'Gap bound packages{", smoothed" if do_smooth  else ""}')
plt.scatter(x, y)
axes.invert_xaxis()

handles, labels = plt.gca().get_legend_handles_labels()

# labels will be the keys of the dict, handles will be values
temp = {k:v for k,v in zip(labels, handles)}

plt.legend(temp.values(), temp.keys(), loc='best', framealpha=1)

plt.show()