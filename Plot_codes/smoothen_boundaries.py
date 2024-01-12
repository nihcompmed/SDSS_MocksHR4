import numpy as np
import helper_functions as hf
from tqdm import tqdm
import pickle
from joblib import Parallel, delayed



sdss_file = '../sdss_all_info.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

tau_u = 22.5

temp = sdss_info.copy()
sdss_save_file = '../sdss_all_info_w_smoothened.p'

for sample in ['sample0', 'sample1', 'sample2']:
    sdss_data = sdss_info['original_undersamples']['percent1'][sample]['cycs']

    temp['original_undersamples']['percent1'][sample]['smoothened_minimal'] = []


    for cyc in tqdm(sdss_data):


        this_boundary = sdss_data[cyc]['minimal_cycs'][0][0]


        this_boundary = hf.local_smoothen_boundary(this_boundary, tau_u)

        temp['original_undersamples']['percent1'][sample]['smoothened_minimal'].append(this_boundary)



pickle.dump(temp, open(sdss_save_file, 'wb'))


mocks_file = '../mocks_all_info.p'
mocks_data = pickle.load(open(mocks_file, 'rb'))
mocks_scale = 0.72

# Convert to Mpc
tau_u = tau_u/0.72

temp = mocks_data.copy()
mocks_save_file = '../mocks_all_info_w_smoothened.p'

mocks_data = mocks_data['mocks_HR4_PH']

def single_mock(mock, cycs):

    #cycs = mocks_data[mock]['cycs']

    smooth_boundaries = []

    for cyc in cycs:
        this_boundary = cycs[cyc]['minimal_cycs'][0]

        this_boundary = hf.local_smoothen_boundary(this_boundary, tau_u)

        smooth_boundaries.append(this_boundary)

    return [mock, smooth_boundaries]


res = Parallel(n_jobs=6, verbose=12)(delayed(single_mock)(mock, mocks_data[mock]['cycs']) for mock in mocks_data.keys())

for r in res:
    mock = r[0]
    smooth_boundaries = r[1]

    temp['mocks_HR4_PH'][mock]['smoothened_minimal'] = smooth_boundaries



pickle.dump(temp, open(mocks_save_file, 'wb'))

