import numpy as np
import helper_functions as hf
from tqdm import tqdm
import pickle
from joblib import Parallel, delayed



sdss_file = '../samethreshold_all_info.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

tau_u = 22.5

temp = sdss_info.copy()
sdss_save_file = '../samethreshold_all_info_w_smoothened.p'

for ss in ['original_undersamples', 'trimmed_undersamples']:

    for percent in ['percent1', 'percent0.95', 'percent0.9', 'percent0.8']:
    
        for sample in ['sample0', 'sample1', 'sample2']:
        
            sdss_data = sdss_info[ss][percent][sample]['cycs']
        
            print(ss, percent, sample)

            temp[ss][percent][sample]['smoothened_minimal'] = []
        
            for cyc in tqdm(sdss_data):
        
                this_boundary = sdss_data[cyc]['minimal_cycs'][0][0]
        
                this_boundary = hf.local_smoothen_boundary(this_boundary, tau_u)
        
                temp[ss][percent][sample]['smoothened_minimal'].append(this_boundary)
    
    
pickle.dump(temp, open(sdss_save_file, 'wb'))



