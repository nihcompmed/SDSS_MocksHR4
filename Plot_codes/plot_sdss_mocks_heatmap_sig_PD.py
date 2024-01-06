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

all_births = []
all_deaths = []

mocks_data = mocks_data['mocks_HR4_PH']

for mock in tqdm(mocks_data):
    pers = mocks_data[mock]['H2_PD']
    
    pers[pers[:,1]==-1, 1] = mock_thresh # mock_thresh is in Mpc

    # Converting all to h^-1 Mpc
    pers = pers*0.72

    diff = pers[:,1] - pers[:,0]

    idxs = np.argwhere(diff >= epsilon).flatten()

    pers = pers[idxs]

    all_births += list(pers[:,0])
    all_deaths += list(pers[:,1])

data = [all_births, all_deaths]

data = np.array(data)

maxx = np.amax(data)
minn = np.amin(data)

kernel = stats.gaussian_kde(data)

#xmin, ymin = np.amin(data, axis=1)
#xmax, ymax = np.amax(data, axis=1)
xmin = minn
ymin = minn
xmax= maxx
ymax = maxx

grid_reso = 0.05

X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
positions = np.vstack([X.ravel(), Y.ravel()])

kernel = stats.gaussian_kde(data)

kernel_estimates = kernel(positions)


Z = np.reshape(kernel_estimates.T, X.shape)

fig, ax = plt.subplots()

plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
          extent=[xmin, xmax, ymin, ymax]\
                  )



# PD of SDSS
sdss_file = '../sdss_all_info.p'
all_info = pickle.load(open(sdss_file, 'rb'))
pers = all_info['original_undersamples']['percent1']['sample1']['H2_PD']

pers[pers[:,1]==-1, 1] = sdss_thresh 

pers[pers[:,1]>sdss_thresh, 1] = sdss_thresh 

diff = pers[:,1] - pers[:,0]

idxs = np.argwhere(diff >= epsilon).flatten()

pers = pers[idxs]

plt.plot(pers[:,0], pers[:,1], marker='x', color='white', markersize=8, lw=0)

plt.xlabel(r'birth ($h^{-1}$ Mpc)', fontsize=16)
plt.ylabel(r'death ($h^{-1}$ Mpc)', fontsize=16)


plt.plot([minn*scale, maxx*scale], [minn*scale, maxx*scale], color='white', ls='--')


plt.savefig('figures/H2_PD_heatmap_SDSS_mocks.jpg', dpi=600)

