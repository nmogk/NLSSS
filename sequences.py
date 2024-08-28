import pickle
from paleobiodb_interface import rv
import requests as req
import json
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime, timezone
from scipy.ndimage import convolve1d

baseUrl = 'https://macrostrat.org/api/'
stagesQuery = f'{baseUrl}defs/intervals?timescale_id=1'

def download_data(opts, header):
    if opts.use_stages:
        res = req.get(stagesQuery)
        data = res.json()['success']['data']
        queries = [x['name'] for x in data]
        x = [x['t_age'] for x in data]
        xt = x
        xb = [x['b_age'] for x in data]
    else:
        step = opts.num
        x = np.arange(0, opts.max_age, step)
        xt = x
        xb = x+step

    y = np.zeros_like(x)
    for i, interval in enumerate(tqdm(x)):
        if opts.use_stages:
            res = req.get(f'{baseUrl}sections?interval_name={queries[i]}{opts.environment_query}')
        else:
            res = req.get(baseUrl+f'sections?age_top={interval}&age_bottom={interval+step}{opts.environment_query}')

        filtered = [x for x in res.json()['success']['data']]
        # Filter out 0 thickness packages
        if opts.filter_zero:
            filtered = [x for x in filtered if x['max_thick'] != '0.00']

        if opts.overlap_type == 'initiate':
            filtered = [x for x in filtered if x['b_age'] <= xb[i] and x['t_age'] < xt[i]]
        elif opts.overlap_type == 'truncate':
            filtered = [x for x in filtered if x['b_age'] > xb[i] and x['t_age'] >= xt[i]]
        elif opts.overlap_type == 'endemic':
            filtered = [x for x in filtered if x['b_age'] <= xb[i] and x['t_age'] >= xt[i]]
        elif opts.overlap_type == 'through':
            filtered = [x for x in filtered if x['b_age'] > xb[i] and x['t_age'] < xt[i]]
        elif opts.overlap_type == 'xupper':
            filtered = [x for x in filtered if x['t_age'] < xt[i]]
        elif opts.overlap_type == 'xlower':
            filtered = [x for x in filtered if x['b_age'] > xb[i]]

        packages = len(filtered) #np.sum([x['col_area'] for x in filtered])
        y[i] = packages
    with open(opts.fname, 'wb') as f:
        dltime = datetime.now(timezone.utc)
        pickle.dump((header, dltime, x, y), f)
    return x, y

def kernel_smooth(y, radius, edge_mode):
    window = 2*radius+1
    kernel = np.array([0.1, .2, 0.4,  .2 ,0.1])#np.ones(window)/window
    return convolve1d(y, kernel, mode=edge_mode)

def plot(x, y, env_title, max_age, do_smooth, use_stages):
    if use_stages:
        sloss_ms = np.array([467, 419, 324, 132, 60])
        peters_ms = np.array([252, 199])
        mogk_ms = np.array([441, 294, 86, 20])
    else:
        sloss_ms = np.array([462, 396, 324, 132, 60])
        peters_ms = np.array([252, 186])
        mogk_ms = np.array([366, 24])

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
    parser.add_argument('-x', action='store_false', dest='compatibility_check', help='Supress data compatibility checking and load the data in the file.')
    parser.add_argument('-t', '--overlap-type', choices=['intersect', 'initiate', 'truncate', 'endemic', 'through', 'xupper', 'xlower'], default='intersect', 
                        help='Types of overlap relationships to count. intersect (default): any part of package overlaps interval; initate: package starts in interval and crosses upper boundary;truncate: package crosses lower boundary and ends in interval; endemic: package is wholly contained in interval; through: package crosses both interval boundaries; xupper: combines initiate and through; xlower: combines truncate and through')

    # Parameters that probably shouldn't be changed, but I don't want to hard code them
    parser.add_argument('--dont-filter-zero', action='store_false', dest='filter_zero', help='Do not remove zero height sediment packages.')
    parser.add_argument('--max-age', type=int, default=540, help='Maximum numeric age bin and horizontal axis value. Default 540Ma')
    parser.add_argument('--kernel-radius', type=int, default=2, metavar='RAD', help='Kernel size will be 2*RAD + 1. Default 2.')
    parser.add_argument('--smoothing-type', choices=['gaussian', 'uniform'], default='gaussian', help='Kernel shape used to perform smoothing.')
    parser.add_argument('--edge-mode', choices=['nearest', 'constant', 'mirror', 'reflect'], default='nearest', help='Edge behavior for smoothing.')
    
    args = parser.parse_args()
    args.use_stages = args.num == 0
    download_settings = dict(bins=args.num, env=args.env, max_age=args.max_age, filter0=args.filter_zero, type=args.overlap_type)

    try:
        with open(args.fname, 'rb') as f:
            header, dltime, x, y = pickle.load(f)

        binfo = str(header["bins"])+" Ma" if header["bins"]!=0 else "stages"

        if args.info:
            print(f'Bins: {binfo}, Env: {header["env"]}, Max Age: {header["max_age"]} Ma, Filter 0 height: {header["filter0"]}, Downloaded: {dltime}')
            return

        if download_settings != header:
            if args.compatibility_check:
                print(f'File exists but contains incompatible data. Bins: {binfo}, Env: {header["env"]}')
                return
            else:
                args.max_age = header["max_age"]
                args.num = header["bins"]
                args.env = header["env"]
                args.use_stages = args.num == 0

    except Exception as e:
        args.environment_query = '' if args.env is None else f'&environ_class={args.env}'
        x, y = download_data(args, download_settings)

    if args.do_smooth:
        y = kernel_smooth(y, args.kernel_radius, args.edge_mode)

    if args.print:
        for a, c in zip((x, y)):
            print(f'{a}, {c}')

    args.env_title = '' if args.env is None else args.env+' '
    plot(x, y, args.env_title, args.max_age, args.do_smooth, args.use_stages)

if __name__ == '__main__':
    main()