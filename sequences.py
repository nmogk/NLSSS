import pickle
from paleobiodb_interface import rv
import requests as req
import json
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

baseUrl = 'https://macrostrat.org/api/'

data_filename = 'packages.pkl'

sloss_ms = np.array([462, 396, 324, 132, 60])
peters_ms = np.array([252, 186])
mogk_ms = np.array([366, 24])

try:
    with open(data_filename, 'rb') as f:
        x, y = pickle.load(f)

except FileNotFoundError:
    step = 6
    x = np.arange(0, 540, step)
    y = np.zeros_like(x)

    with open(data_filename, 'wb') as f:
        for i, interval in enumerate(tqdm(x)):
            # res = req.get(baseUrl+f'sections?interval_name={interval[rv.NAME]}&environ_class=marine')
            res = req.get(baseUrl+f'sections?age_top={interval}&age_bottom={interval+step}&environ_class=marine')

            # Filter out 0 thickness packages
            filtered = [x for x in res.json()['success']['data'] if x['max_thick'] != '0.00']
            packages = len(filtered)#np.sum([x['col_area'] for x in filtered]) #len(filtered)
            y[i] = packages
        pickle.dump((x,y), f)

window = 5
kernel = np.array([0.1, .2, 0.4,  .2 ,0.1])#np.ones(window)/window
y = np.convolve(y, kernel, mode='same')
# print(x)
# print(y)

plt.title('Number of gap bound packages')
for bdry in sloss_ms:
    plt.axvline(bdry, color='gray', linestyle='--')
for bdry in peters_ms:
    plt.axvline(bdry, color='blue', linestyle='--')
for bdry in mogk_ms:
    plt.axvline(bdry, color='gold', linestyle=':')
plt.plot(x,y)
plt.scatter(x, y)
plt.gca().invert_xaxis()
plt.show()