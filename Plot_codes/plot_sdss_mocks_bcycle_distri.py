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

# The computation for SDSS were done till 36 h^-1 Mpc in an attempt to get persistence of all
# features
# But it is limited to 30 here for plotting purpose to match mocks threshold of 30 h^-1 Mpc
sdss_thresh = 30 # In h^-1 Mpc

# Scale to convert from Mpc to h^-1 Mpc
mocks_scale = 0.72
mock_thresh = round(30/0.72,2) # In Mpc

mocks_file = '../mocks_all_info.p'
mocks_data = pickle.load(open(mocks_file, 'rb'))

n_cyc = []

mocks_data = mocks_data['mocks_HR4_PH']

for mock in tqdm(mocks_data):
    cycles = mocks_data[mock]['cycs']
    
    n_cyc.append(len(cycles))


sns.kdeplot(n_cyc)


# PD of SDSS
sdss_file = '../sdss_all_info.p'
all_info = pickle.load(open(sdss_file, 'rb'))
cycs = all_info['original_undersamples']['percent1']['sample1']['cycs']

sdss_n_sig = len(cycs)

plt.gca().axvline(x=sdss_n_sig, color='red', ls='--', lw=2, label=r'\#tight rep. in SDSS')

plt.xlabel(r'Number of tight representatives computed for mock galaxies', fontsize=16)
plt.ylabel('Density', fontsize=16)

plt.legend(fontsize=14)


plt.savefig('figures/n_bcycles_korea_mocks_sdss_H2.jpg', dpi=600)





