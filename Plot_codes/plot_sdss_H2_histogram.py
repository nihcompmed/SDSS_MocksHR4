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
sdss_thresh = 36 # In h^-1 Mpc

# PD of SDSS
sdss_file = '../sdss_all_info.p'
all_info = pickle.load(open(sdss_file, 'rb'))

pers = all_info['original_undersamples']['percent1']['sample1']['H2_PD']

pers[pers[:,1]==-1, 1] = sdss_thresh 

diff = pers[:,1] - pers[:,0]

plt.hist(diff, bins=75)

plt.axvline(x=epsilon, color='darkorange', ls='--', lw=2, label='persistence threshold')

plt.yscale('log')
plt.xlabel(r'persistence ($h^{-1}$ Mpc)', fontsize=16)
plt.ylabel(r'counts', fontsize=16)
plt.legend(fontsize=14)

plt.savefig('figures/SDSS_H2_pers_distri.jpg', dpi=600)

