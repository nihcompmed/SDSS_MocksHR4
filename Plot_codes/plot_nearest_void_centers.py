import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from tqdm import tqdm
import helper_functions as hf
from sklearn.neighbors import KernelDensity

# Setting for kde
X_plot = np.linspace(0, 150, 150)[:, np.newaxis]
bandwidth = 2

sdss_file = '../sdss_all_info.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

sdss_samples_kde = []


for sample in ['sample0', 'sample1', 'sample2']:
    sdss_data = sdss_info['original_undersamples']['percent1'][sample]['cycs']
    locs = sdss_info['original_undersamples']['percent1'][sample]['locs']

    boundaries = []
    
    for cyc in tqdm(sdss_data):
        this_boundary = sdss_data[cyc]['minimal_cycs'][0][0]
        boundaries.append(this_boundary)
    
    dists = hf.get_nearest_void_center_distances(boundaries)

    plt.hist(dists, bins=20, alpha=0.5)
    plt.xlabel(r'distance to nearest void center ($h^{-1}$ Mpc)', fontsize=14)
    plt.ylabel(r'counts', fontsize=16, color='tab:blue')
    plt.gca().tick_params(axis='y', labelcolor='tab:blue')

    ax2 = plt.gca().twinx()

    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(dists[:, np.newaxis])
    log_dens = kde.score_samples(X_plot)

    dens = np.exp(log_dens)

    sdss_samples_kde.append(dens)

    ax2.plot(X_plot, dens, color='red', zorder=1000)
    ax2.set_ylabel('Density', color='red', fontsize=16)
    ax2.tick_params(axis='y', labelcolor='red')

    plt.tight_layout()

    plt.savefig(f'figures/SDSS_{sample}_nearest_void_centers.jpg', dpi=600)
    plt.cla()
    plt.clf()

plt.cla()
plt.clf()

# Sample 1???
plt.plot(X_plot, sdss_samples_kde[1], color='r', alpha=0.6, zorder=4000)


mocks_file = '../mocks_all_info.p'

mocks_data = pickle.load(open(mocks_file, 'rb'))

mocks_data = mocks_data['mocks_HR4_PH']
# To convert to h^-1 Mpc
mocks_scale = 0.72


all_mocks_kde = []

for mock in tqdm(mocks_data):

    cycs = mocks_data[mock]['cycs']

    boundaries = []

    for cyc in cycs:
        # convert to h^-1 Mpc
        this_boundary = mocks_data[mock]['cycs'][cyc]['minimal_cycs'][0]*mocks_scale
        boundaries.append(this_boundary)

    dists = hf.get_nearest_void_center_distances(boundaries)

    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(dists[:, np.newaxis])
    log_dens = kde.score_samples(X_plot)

    dens = np.exp(log_dens)

    plt.plot(X_plot, dens, color='turquoise', alpha=0.4)

    all_mocks_kde.append(dens)

all_mocks_kde = np.array(all_mocks_kde)

ptiles = np.percentile(all_mocks_kde, [10, 25, 50, 75, 90], axis=0)

## 25 ptile
#ptile25 = np.vstack((ptiles[1],ptiles[-2]))
#
## 10 ptile
#ptile10 = np.vstack((ptiles[0],ptiles[-1]))

# plot median
plt.plot(np.squeeze(X_plot), ptiles[2], color='blue', zorder=4000, label='mocks median')

plt.fill_between(np.squeeze(X_plot), ptiles[1], ptiles[-2]\
                            , zorder=2000, alpha=0.4, color='blue', label=r'mocks 25 to 75%-tile')
plt.fill_between(np.squeeze(X_plot), ptiles[0], ptiles[-1]\
                            , zorder=2000, alpha=0.4, color='orange', label=r'mocks 10 to 90%-tile')


old_handles, labels = plt.gca().get_legend_handles_labels()
legend_elements = [Line2D([0], [0], color='r', lw=2, label='SDSS') ]

plt.legend(handles=old_handles + legend_elements)


plt.xlabel(r'distance to nearest void center ($h^{-1}$ Mpc)', fontsize=14)
plt.ylabel('Density', fontsize=16)

plt.savefig('figures/SDSS_mocks_void_centers_kde.jpg', dpi=600)










#    all_gor = []
#    for r_vals, this_gor in this_sample:
#        all_gor.append(this_gor)
#        if sample == 'sample0':
#            plt.plot(r_vals, this_gor, alpha=0.2, color='red')
#
#    all_gor = np.array(all_gor)
#
#    gor_med = np.median(all_gor, axis=0)
#
#    sdss_sample_medians.append(gor_med)
#
#    if sample == 'sample0':
#        plt.plot(r_vals, gor_med, color='red', lw=2, label='median')
#        plt.legend()
#
#        plt.xlabel(r'$r$ ($h^{-1}$ Mpc)', fontsize=16)
#        plt.ylabel(r'$g(r)$', fontsize=16)
#
#        plt.savefig('figures/SDSS_sample0_gor.jpg', dpi=600)
#
#plt.cla()
#plt.clf()
