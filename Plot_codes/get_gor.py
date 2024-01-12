import numpy as np
import helper_functions as hf
from tqdm import tqdm
import pickle



sdss_file = '../sdss_all_info_w_smoothened.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

sdss_gor = dict()

for sample in ['sample0', 'sample1', 'sample2']:
    sdss_data = sdss_info['original_undersamples']['percent1'][sample]['smoothened_minimal']
    locs = sdss_info['original_undersamples']['percent1'][sample]['locs']

    sdss_gor[sample] = []
    
    for this_boundary in tqdm(sdss_data):
        r_vals, gor = hf.get_gor(this_boundary, locs)
        sdss_gor[sample].append((r_vals, gor))

pickle.dump(sdss_gor, open('sdss_gor.p', 'wb'))

mocks_gor = dict()

mocks_file = '../mocks_all_info_w_smoothened.p'
mocks_data = pickle.load(open(mocks_file, 'rb'))
mocks_data = mocks_data['mocks_HR4_PH']
mocks_scale = 0.72

for mock in tqdm(mocks_data):

    locs = mocks_data[mock]['locs']*mocks_scale

    mocks_gor[mock] = []

    for this_boundary in mocks_data[mock]['smoothened_minimal']:
        # convert to h^-1 Mpc
        this_boundary = this_boundary*mocks_scale
        r_vals, gor = hf.get_gor(this_boundary, locs)
        mocks_gor[mock].append((r_vals, gor))


pickle.dump(mocks_gor, open('mocks_gor.p', 'wb'))




