import numpy as np
import helper_functions as hf
from tqdm import tqdm
import pickle



sdss_file = '../sdss_all_info.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

sdss_gor = dict()

for sample in ['sample0', 'sample1', 'sample2']:
    sdss_data = sdss_info['original_undersamples']['percent1'][sample]['cycs']
    locs = sdss_info['original_undersamples']['percent1'][sample]['locs']

    sdss_gor[sample] = []
    
    for cyc in tqdm(sdss_data):
        this_boundary = sdss_data[cyc]['minimal_cycs'][0][0]
        r_vals, gor = hf.get_gor(this_boundary, locs)
        sdss_gor[sample].append((r_vals, gor))

pickle.dump(sdss_gor, open('sdss_gor.p', 'wb'))

mocks_gor = dict()

mocks_file = '../mocks_all_info.p'
mocks_data = pickle.load(open(mocks_file, 'rb'))
mocks_data = mocks_data['mocks_HR4_PH']
mocks_scale = 0.72

for mock in tqdm(mocks_data):

    cycs = mocks_data[mock]['cycs']

    locs = mocks_data[mock]['locs']*mocks_scale

    mocks_gor[mock] = []

    for cyc in cycs:
        # convert to h^-1 Mpc
        this_boundary = mocks_data[mock]['cycs'][cyc]['minimal_cycs'][0]*mocks_scale
        r_vals, gor = hf.get_gor(this_boundary, locs)
        mocks_gor[mock].append((r_vals, gor))


pickle.dump(mocks_gor, open('mocks_gor.p', 'wb'))




