import pickle
from paleobiodb_interface import rv
import requests as req
import json
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime, timezone

baseUrl = 'https://macrostrat.org/api/'
stagesQuery = f'{baseUrl}defs/intervals?timescale_id=1'

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

# try:
#     with open(data_filename, 'rb') as f:
#         x, y = pickle.load(f)

# except:
#     pass

def download_data(opts, header):
    if opts.use_stages:
        res = req.get(stagesQuery)
        data = res.json()['success']['data']
        queries = [x['name'] for x in data]
        x = [x['t_age'] for x in data]
    else:
        step = opts.num
        x = np.arange(0, opts.max_age, step)

    y = np.zeros_like(x)
    with open(opts.fname, 'wb') as f:
        dltime = datetime.now(timezone.utc)
        for i, interval in enumerate(tqdm(x)):
            if opts.use_stages:
                res = req.get(f'{baseUrl}sections?interval_name={queries[i]}{opts.environment_query}')
            else:
                res = req.get(baseUrl+f'sections?age_top={interval}&age_bottom={interval+step}{opts.environment_query}')

            # Filter out 0 thickness packages
            if opts.dont_filter_zero:
                filtered = [x for x in res.json()['success']['data']]
            else:
                filtered = [x for x in res.json()['success']['data'] if x['max_thick'] != '0.00']
            packages = len(filtered) #np.sum([x['col_area'] for x in filtered]) #len(filtered)
            y[i] = packages
        pickle.dump((header, dltime, x, y), f)
    return x, y

def kernel_smooth(y, radius, edge_mode):
    window = 2*radius+1
    kernel = np.array([0.1, .2, 0.4,  .2 ,0.1])#np.ones(window)/window
    return np.convolve(y, kernel, mode=edge_mode)

def plot(x, y, env_title, max_age, do_smooth, use_stages):
    plt.title(f'Number of {env_title}gap bound packages{", smoothed" if do_smooth  else ""}')
    axes = plt.gca()
    axes.set_xlim([0, max_age])
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

def main():
    import argparse

    parser = argparse.ArgumentParser(description= 'Find megasequences using macrostrat.')

    # Core options
    parser.add_argument('fname', metavar='FILENAME', help='Name of data file to analyze or cache downloads.')
    parser.add_argument('-n', '--num', nargs='?', type=int, default=0, const=6, metavar='BIN', 
                        help='Use numeric age bins instead of stages. Optionally set BIN size, default 6Ma.')
    parser.add_argument('-e', '--env', choices=['marine', 'non-marine'], help='Only use packages from specified environment.')
    parser.add_argument('-s', '--smooth', action='store_true', dest='do_smooth', help='Apply a windowed filter to the output.')
    parser.add_argument('-p', '--print', action='store_true', help='Print results.')
    parser.add_argument('-i', '--info', action='store_true', help='Summarize the file download header and exit.')

    # Parameters that probably shouldn't be changed, but I don't want to hard code them
    parser.add_argument('--dont-filter-zero', action='store_true', help='Do not remove zero height sediment packages.')
    parser.add_argument('--max-age', type=int, default=540, help='Maximum numeric age bin and horizontal axis value. Default 540Ma')
    parser.add_argument('--kernel-radius', type=int, default=2, metavar='RAD', help='Kernel size will be 2*RAD + 1. Default 2.')
    parser.add_argument('--smoothing-type', choices=['gaussian', 'uniform'], default='gaussian', help='Kernel shape used to perform smoothing.')
    parser.add_argument('--edge-mode', choices=['same', 'nearest'], default='same', help='Edge behavior for smoothing.')
    
    args = parser.parse_args()
    args.use_stages = args.num == 0
    args.env_title = '' if args.env is None else args.env+' '
    args.environment_query = '' if args.env is None else f'&environ_class={args.env}'

    download_settings = dict(bins=args.num, env=args.env, max_age=args.max_age, filter0=not args.dont_filter_zero)

    try:
        with open(args.fname, 'rb') as f:
            header, dltime, x, y = pickle.load(f)

        binfo = str(header["bins"])+" Ma" if header["bins"]!=0 else "stages"

        if args.info:
            print(f'Bins: {binfo}, Env: {header["env"]}, Max Age: {header["max_age"]} Ma, Filter 0 height: {header["filter0"]}, Downloaded: {dltime}')
            return

        if download_settings != header:
            print(f'File exists but contains incompatible data. Bins: {binfo}, Env: {header["env"]}')
            return

    except Exception as e:
        print(e)
        x, y = download_data(args, download_settings)

    if args.do_smooth:
        y = kernel_smooth(y, args.kernel_radius, args.edge_mode)

    if args.print:
        for a, c in zip((x, y)):
            print(f'{a}, {c}')

    plot(x, y, args.env_title, args.max_age, args.do_smooth, args.use_stages)

if __name__ == '__main__':
    main()