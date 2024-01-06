import numpy as np
import matplotlib.pyplot as plt
import os
import math
import seaborn as sns
from scipy import stats
import pickle
from tqdm import tqdm
plt.rcParams['text.usetex'] = True


# Units in h^-1 Mpc
tau_u = 22.5
epsilon = 7.5
scale = 1

# This is strictly for plotting purpose to match mocks threshold of 30 h^-1 Mpc
sdss_thresh = 30 # In h^-1 Mpc

# Scale to convert from Mpc to h^-1 Mpc
mocks_scale = 0.72
mock_thresh = round(30/0.72,2) # In Mpc

mocks_file = '../mocks_all_info.p'
mocks_data = pickle.load(open(mocks_file, 'rb'))

n_sig = []

mocks_data = mocks_data['mocks_HR4_PH']

for mock in tqdm(mocks_data):
    pers = mocks_data[mock]['H2_PD']
    
    pers[pers[:,1]==-1, 1] = mock_thresh # mock_thresh is in Mpc

    # Converting all to h^-1 Mpc
    pers = pers*0.72

    diff = pers[:,1] - pers[:,0]

    idxs = np.argwhere(diff >= epsilon).flatten()

    n_sig.append(len(idxs))


sns.kdeplot(n_sig)


# PD of SDSS
sdss_file = '../sdss_all_info.p'
all_info = pickle.load(open(sdss_file, 'rb'))
pers = all_info['original_undersamples']['percent1']['sample1']['H2_PD']

pers[pers[:,1]==-1, 1] = sdss_thresh 

pers[pers[:,1]>sdss_thresh, 1] = sdss_thresh 

diff = pers[:,1] - pers[:,0]

idxs = np.argwhere(diff >= epsilon).flatten()

sdss_n_sig = len(idxs)

plt.gca().axvline(x=sdss_n_sig, color='red', ls='--', lw=2, label=r'\#sig. in SDSS')

plt.xlabel(r'Number of H$_2$ significant features in mock galaxies', fontsize=16)
plt.ylabel('Density', fontsize=16)

plt.legend(fontsize=14)


plt.savefig('figures/n_sig_korea_mocks_sdss_H2.jpg', dpi=600)


