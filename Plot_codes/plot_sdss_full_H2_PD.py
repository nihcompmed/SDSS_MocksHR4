import numpy as np
import matplotlib.pyplot as plt
import os
import math
import seaborn as sns
from scipy import stats
import pickle



def plot_PD(pers, epsilon, tau_u, filename=None, scale=1):


    epsilon = epsilon*scale
    tau_u = tau_u*scale

    undead = np.argwhere(pers[:,1] == -1).flatten()
    dead = np.argwhere(pers[:,1] != -1).flatten()

    pers_undead = pers[undead]
    pers_dead = pers[dead]

    if pers_undead.ndim == 1:
        pers_undead = np.reshape(pers_undead, (1, 2))
    if pers_dead.ndim == 1:
        pers_dead = np.reshape(pers_dead, (1, 2))


    

    pers_dead = pers_dead*scale
    diff_dead = pers_dead[:,1] - pers_dead[:,0]
    plt.scatter(pers_dead[:,0], pers_dead[:,1], color='blue', alpha=0.4, s=15)
    #plt.scatter(pers_dead[:,0], diff_dead, color='blue', alpha=0.4)

    minn = np.amin(pers_dead)
    maxx = np.amax(pers_dead)

    undead_line = 1*np.amax(pers)
    pers_undead[:,1] = undead_line

    pers_undead = pers_undead*scale
    diff_undead = pers_undead[:,1] - pers_undead[:,0]
    plt.scatter(pers_undead[:,0], pers_undead[:,1], color='red', alpha=0.4, s=15)
    #plt.scatter(pers_undead[:,0], diff_undead, color='red', alpha=0.4)
    plt.gca().axhline(y=undead_line*scale, alpha=0.5, ls='--')

    #locs, labels = plt.yticks() 
    #locs = list(locs)
    #locs.append(undead_line*scale)
    #labels += [r'$\infty$']
    #plt.yticks(ticks=locs, labels=labels)


    plt.plot([minn, maxx-epsilon], [minn+epsilon, maxx], ls='--', color='orange', lw=2, label='persistence threshold')

    plt.gca().axvline(x=tau_u, color='black', ls='--', lw=2, label='birth threshold')
    #plt.gca().axhline(y=epsilon, color='orange', ls='--', lw=2, label='persistence threshold')


    plt.plot([minn, maxx], [minn, maxx], ls='--', color='black')

    plt.xlabel(r'birth $(h^{-1} \, {\rm Mpc})$', fontsize=16)
    plt.ylabel(r'death $(h^{-1} \, {\rm Mpc})$', fontsize=16)

    plt.legend()
    
    #plt.savefig('figures/'+galaxy_id+'_H2_PD.png')
    #plt.show()
    #exit()

    #plt.cla()
    #plt.clf()
    if filename:
        plt.savefig(filename, dpi=600)


    plt.cla()
    plt.clf()



# Units in h^-1 Mpc
tau_u = 22.5
epsilon = 7.5
scale = 1


# Plot PD of SDSS
sdss_file = '../sdss_all_info.p'
all_info = pickle.load(open(sdss_file, 'rb'))

sdss_data = all_info['original_undersamples']['percent1']['sample1']['H2_PD']


plot_PD(sdss_data, epsilon, tau_u, 'figures/SDSS_sample1_H2_PD.jpg', scale)



##this_n_sig = plot_PD(korea_sdss_file, epsilon, tau_u, 'SDSS_H2.pdf')
#
#sdss_pers = np.loadtxt(korea_sdss_file, delimiter=',')
#sdss_pers[sdss_pers[:,1]==-1, 1] = thresh
#diff = sdss_pers[:,1] - sdss_pers[:,0]
#idxs = np.argwhere(diff >= epsilon).flatten()
#sdss_pers = sdss_pers[idxs]
#
#
#
#
#ff = glob.glob(dirr+'/*')
#
#
#all_birth = []
#all_death = []
#
#count = 0
#for dirr in ff:
#
#    # ADD TO IGNORE SDSS
#
#    pers_file = dirr + '/H2_pers_data.txt'
#
#    print(count, end='\r')
#    count += 1
#    #if count == 10:
#    #    break
#
#    if not os.path.isfile(pers_file):
#        continue
#
#    pers = np.loadtxt(pers_file, delimiter=',')
#
#    pers = pers[pers[:,0] <= tau_u]
#
#    pers[pers[:,1]==-1, 1] = thresh
#
#    diff = pers[:,1] - pers[:,0]
#
#    idxs = np.argwhere(diff >= epsilon).flatten()
#
#    pers = pers[idxs]
#
#    all_birth += list(pers[:,0])
#    all_death += list(pers[:,1])
#
#
#data = [all_birth, all_death]
#
#data = np.array(data)
#
#maxx = np.amax(data)
#minn = np.amin(data)
#
#kernel = stats.gaussian_kde(data)
#
#
##xmin, ymin = np.amin(data, axis=1)
##xmax, ymax = np.amax(data, axis=1)
#xmin = minn
#ymin = minn
#xmax= maxx
#ymax = maxx
#
#grid_reso = 0.00005
#
#X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
#positions = np.vstack([X.ravel(), Y.ravel()])
#
#kernel = stats.gaussian_kde(data)
#
#kernel_estimates = kernel(positions)
#
#
#Z = np.reshape(kernel_estimates.T, X.shape)
#
#fig, ax = plt.subplots()
#
#plt.imshow(np.rot90(Z*scale), cmap=plt.cm.gnuplot,
#          extent=[xmin*scale, xmax*scale, ymin*scale, ymax*scale]\
#                  )
#
#plt.plot(sdss_pers[:,0]*scale, sdss_pers[:,1]*scale, marker='x', color='white', markersize=8, lw=0)
#
#plt.xlim([xmin*scale, xmax*scale])
#plt.ylim([ymin*scale, ymax*scale])
#plt.xlabel('birth', fontsize=16)
#plt.ylabel('death', fontsize=16)
#
#
#plt.plot([minn*scale, maxx*scale], [minn*scale, maxx*scale], color='white', ls='--')
#
#
#plt.savefig('H2_PD_heatmatp_SDSS_mocks.pdf', dpi=600)
#
#
#
#exit()
#
##kernel = stats.gaussian_kde(data)
##
##X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
##positions = np.vstack([X.ravel(), Y.ravel()])
#
#
#
#
#
##sns.violinplot(all_n_sig, color='turquoise')
#
#
#
#
#sdss_n_sig = get_n_sig(korea_sdss_file, epsilon, tau_u)
#
#plt.gca().axhline(y=sdss_n_sig, color='red', ls='--', lw=2, label='#sig. in SDSS')
#
#plt.ylabel(r'Number of H$_2$ significant features in mock galaxies')
#
#plt.legend()
#
#plt.savefig('n_sig_korea_mocks_sdss_H2.pdf', dpi=600)
#
#plt.cla()
#plt.clf()
#
#all_n_sig = []
#
#count = 0
#for dirr in ff:
#
#    pers_file = dirr + '/H1_pers_data.txt'
#
#    if not os.path.isfile(pers_file):
#        continue
#
#    galaxy_id = dirr.split('/')[1]
#
#    this_n_sig = get_n_sig(pers_file, epsilon, tau_u)
#
#    all_n_sig.append(this_n_sig)
#
#sns.violinplot(all_n_sig, color='turquoise')
#
#korea_sdss_file = 'results/Korea_SDSS_amd128/thresh0.01_b0.0075_H1_pers_data.txt'
#
#sdss_n_sig = get_n_sig(korea_sdss_file, epsilon, tau_u)
#
#plt.gca().axhline(y=sdss_n_sig, color='red', ls='--', lw=2, label='#sig. in SDSS')
#
#plt.ylabel(r'Number of H$_1$ significant features in mock galaxies')
#
#plt.legend()
#
#plt.savefig('n_sig_korea_mocks_sdss_H1.pdf', dpi=600)
#
#
#
#
#
#
#
#
#
#

