import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import math
import seaborn as sns
from scipy import stats
import pickle
import helper_functions as hf
from tqdm import tqdm
from sklearn.neighbors import KernelDensity
plt.rcParams['text.usetex'] = True

# Setting for kde
X_plot = np.linspace(0, 80, 100)[:, np.newaxis]
bandwidth = 2


sdss_file = '../sdss_all_info_w_smoothened.p'
all_info = pickle.load(open(sdss_file, 'rb'))

all_sdss_minn = []
all_sdss_med = []
all_sdss_maxx = []

for sample in ['sample0', 'sample1', 'sample2']:

    # Only plot sample1 for clean plot
    if sample != 'sample1':
        continue

    sdss_data = all_info['original_undersamples']['percent1'][sample]['smoothened_minimal']
    
    sdss_reff = []
    
    for this_boundary in sdss_data:
        this_reff = hf.get_reff_convexhull(this_boundary)
        #this_reff = hf.get_reff_delaunay(this_boundary)
        sdss_reff.append(this_reff)

    sdss_reff = np.array(sdss_reff)

    all_sdss_minn.append(np.amin(sdss_reff))
    all_sdss_maxx.append(np.amax(sdss_reff))
    all_sdss_med.append(np.median(sdss_reff))
    
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(sdss_reff[:, np.newaxis])
    log_dens = kde.score_samples(X_plot)
    plt.plot(X_plot, np.exp(log_dens), color='red', zorder=1000)

    #sns.kdeplot(sdss_reff, bw_method=1, zorder=1000)

#plt.show()
#exit()

mocks_file = '../mocks_all_info_w_smoothened.p'

mocks_data = pickle.load(open(mocks_file, 'rb'))

mocks_data = mocks_data['mocks_HR4_PH']
# To convert to h^-1 Mpc
mocks_scale = 0.72

skip = 0

all_mocks = []

all_mocks_minn = []
all_mocks_maxx = []
all_mocks_med = []

count = 0

for mock in tqdm(mocks_data):

    ##########
    ## Testing
    #if count == 10:
    #    break
    #count += 1
    ##########

    cycs = mocks_data[mock]['cycs']
    flag = 1
    for cyc in cycs:
        n_sig = mocks_data[mock]['cycs'][cyc]['n_sig']
        if n_sig != 1:
            flag = 0
    if flag:
        mock_reff = []
        for this_boundary in mocks_data[mock]['smoothened_minimal']:
            # convert to h^-1 Mpc
            this_boundary = this_boundary*mocks_scale
            this_reff = hf.get_reff_convexhull(this_boundary)
            mock_reff.append(this_reff)

        mock_reff = np.array(mock_reff)
        all_mocks_minn.append(np.amin(mock_reff))
        all_mocks_maxx.append(np.amax(mock_reff))
        all_mocks_med.append(np.median(mock_reff))

        #sns.kdeplot(mock_reff, bw_method=1, alpha=0.4, color='turquoise')
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(mock_reff[:, np.newaxis])
        log_dens = kde.score_samples(X_plot)
        dens = np.exp(log_dens)
        all_mocks.append(dens)
        plt.plot(X_plot, dens, color='turquoise', alpha=0.2)
    else:
        skip += 1

print(f'Skipping {skip} mocks that did not have exactly one feature in the computed tight representatives')

all_mocks = np.array(all_mocks)



mocks_med = np.median(all_mocks, axis=0)
ptiles = np.percentile(all_mocks, [10, 25, 50, 75, 90], axis=0)

## 25 ptile
#ptile25 = np.vstack((ptiles[1],ptiles[-2]))
#
## 10 ptile
#ptile10 = np.vstack((ptiles[0],ptiles[-1]))

# plot median
plt.plot(np.squeeze(X_plot), ptiles[2], color='blue', zorder=2000, label='mocks median')

plt.fill_between(np.squeeze(X_plot), ptiles[1], ptiles[-2]\
                            , zorder=2000, alpha=0.4, color='blue', label=r'mocks 25 to 75%-tile')
plt.fill_between(np.squeeze(X_plot), ptiles[0], ptiles[-1]\
                            , zorder=2000, alpha=0.4, color='orange', label=r'mocks 10 to 90%-tile')


old_handles, labels = plt.gca().get_legend_handles_labels()
legend_elements = [Line2D([0], [0], color='r', lw=2, label='SDSS') ]

plt.legend(handles=old_handles + legend_elements)


plt.xlabel(r'$R_{\rm eff}$ ($h^{-1}$ Mpc)', fontsize=18)
plt.ylabel('Density', fontsize=18)

    
plt.savefig('figures/SDSS_mocks_Reff.jpg', dpi=600)

plt.cla()
plt.clf()


plt.scatter([0]*len(all_mocks_minn), all_mocks_minn, alpha=0.4, color='turquoise')
plt.scatter([1]*len(all_mocks_med), all_mocks_med, alpha=0.4, color='turquoise')
plt.scatter([2]*len(all_mocks_maxx), all_mocks_maxx, alpha=0.4, color='turquoise')

plt.scatter([0]*len(all_sdss_minn), all_sdss_minn, color='red', alpha=0.6)
plt.scatter([1]*len(all_sdss_med), all_sdss_med, color='red', alpha=0.6)
plt.scatter([2]*len(all_sdss_maxx), all_sdss_maxx, color='red', alpha=0.6)

plt.ylabel(r'$R_{\rm eff}$ ($h^{-1}$ Mpc)', fontsize=18)

plt.xlim([-0.75, 2.75])

plt.xticks([0, 1, 2], ['min', 'median', 'max'], fontsize=18)


plt.savefig('figures/SDSS_mocks_reff_minmaxmed.jpg', dpi=600)


## Scale to convert from Mpc to h^-1 Mpc
#mocks_scale = 0.72
#mock_thresh = round(30/0.72,2) # In Mpc
#
#mocks_file = '../mocks_all_info.p'
#mocks_data = pickle.load(open(mocks_file, 'rb'))
#
#n_cyc = []
#
#mocks_data = mocks_data['mocks_HR4_PH']
#
#for mock in tqdm(mocks_data):
#    cycles = mocks_data[mock]['cycs']
#    
#    n_cyc.append(len(cycles))
#
#
#sns.kdeplot(n_cyc)
#
#
## PD of SDSS
#sdss_file = '../sdss_all_info.p'
#all_info = pickle.load(open(sdss_file, 'rb'))
#cycs = all_info['original_undersamples']['percent1']['sample0']['cycs']
#
#sdss_n_sig = len(cycs)
#
#plt.gca().axvline(x=sdss_n_sig, color='red', ls='--', lw=2, label=r'\#tight rep. in SDSS')
#
#plt.xlabel(r'Number of tight representatives computed for mock galaxies', fontsize=16)
#plt.ylabel('Density', fontsize=16)
#
#plt.legend(fontsize=14)
#
#
#plt.savefig('figures/n_bcycles_korea_mocks_sdss_H2.pdf', dpi=600)
#


