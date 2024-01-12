import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import math
import seaborn as sns
import pickle
import helper_functions as hf
from tqdm import tqdm
import itertools as it
from scipy.stats import mannwhitneyu


sdss_file = '../samethreshold_all_info_w_smoothened.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

locs = sdss_info['original_undersamples']['percent1']['sample0']['locs']

# Plot #sig features

colors = dict()
markers = dict()

all_centers = dict()

markers_list = ['o', 'x']
colors_list = ['#D81B60', '#1E88E5', '#FFC107']

idx = 0
for typ in sdss_info:
    markers[typ] = markers_list[idx]
    idx2 = 0
    for sample in sdss_info[typ]['percent1']:
        colors[sample] = colors_list[idx2]
        idx2 += 1
    
    idx += 1

all_typ = list(sdss_info.keys())
all_per = list(sdss_info[all_typ[0]].keys())
all_sample = list(sdss_info[all_typ[0]][all_per[0]].keys())

row = 0
col = 0

# Parameters in h^-1 Mpc
tau_u = 22.5
epsilon = 7.5
thresh = 35

percents = [1, 0.95, 0.9, 0.8]

maxx = 0

plot_PD_typ = []
plot_PD_per = []
plot_PD_val = []

plot_tight_typ = []
plot_tight_per = []
plot_tight_val = []


for typ in sdss_info:

    this_typ = sdss_info[typ]

    for per in all_per:

        if per == 'percent1':
            percent = 100
        elif per == 'percent0.95':
            percent = 95
        elif per == 'percent0.9':
            percent = 90
        elif per == 'percent0.8':
            percent = 80

        print(per, percent)


        this_tau_u = round(tau_u, 2)
        this_epsilon = round(epsilon, 2)
        this_thresh = round(thresh, 2)

        this_per = this_typ[per]

        idx = 0

        for sample in this_per:

            color = colors_list[idx]

            this_sample = this_per[sample]

            this_sample_H2 = this_sample['H2_PD']

            this_sample_H2 = this_sample_H2[this_sample_H2[:,0] <= this_tau_u]

            this_sample_H2[this_sample_H2[:,1]==-1, 1] = this_thresh

            this_sample_H2 = this_sample_H2[this_sample_H2[:,1] - this_sample_H2[:,0] >= this_epsilon]

            # Assuming more than 1
            this_n_sig_PD = this_sample_H2.shape[0]

            this_n_sig_found = len(this_sample['smoothened_minimal'])

            if typ == 'original_undersamples':
                plot_PD_typ.append('original SDSS')
                plot_tight_typ.append('original SDSS')
            elif typ == 'trimmed_undersamples':
                plot_PD_typ.append('trimmed SDSS')
                plot_tight_typ.append('trimmed SDSS')

            plot_PD_per.append(str(percent))
            plot_PD_val.append(this_n_sig_PD)

            plot_tight_per.append(str(percent))
            plot_tight_val.append(this_n_sig_found)


            #found_pers = np.array(found_pers)
            #print(found_pers)
            #axs[row, col].scatter(found_pers[:,0], found_pers[:,1]\
            #                    , marker='x', color = color, alpha=0.5)
            
            idx += 1
                    
fig, axs = plt.subplots(1, 2, sharey=True)
sns.swarmplot(x=plot_PD_per, y=plot_PD_val, hue=plot_PD_typ, size=10, ax=axs[0])

sns.swarmplot(x=plot_tight_per, y=plot_tight_val, hue=plot_tight_typ, size=10, ax=axs[1])

axs[0].set_ylabel('Number of significant features', fontsize=12)

axs[0].set_xlabel('Random sampling percentage', fontsize=12)
axs[1].set_xlabel('Random sampling percentage', fontsize=12)

axs[0].set_title('PD computation', fontsize=12)
axs[1].set_title('Tight representatives computation', fontsize=12)

plt.savefig('figures/num_sig_PD_tight_samethresh.jpg', dpi=600)

plt.cla()
plt.clf()
plt.close()

# Compare centroids (plot xyz)

minns = np.min(locs, axis=0)
maxxs = np.max(locs, axis=0)


colors = dict()
markers = dict()

all_centers = dict()

markers_list = ['o', 'x']

markers_dict = {'original_undersamples':'o', 'trimmed_undersamples':'x'}
colors_list = ['#D81B60', '#1E88E5', '#FFC107']

idx = 0
for typ in sdss_info:
    markers[typ] = markers_list[idx]
    idx2 = 0
    for sample in sdss_info[typ]['percent1']:
        colors[sample] = colors_list[idx2]
        idx2 += 1
    
    idx += 1

all_typ = list(sdss_info.keys())
all_per = list(sdss_info[all_typ[0]].keys())
all_sample = list(sdss_info[all_typ[0]][all_per[0]].keys())

#fig, axs = plt.subplots(4, 2, figsize=(10,8), sharex=True, sharey=True)

row = 0
col = 0

percents = [1, 0.95, 0.9, 0.8]

fig, axs = plt.subplots(4, 3, figsize=(10,8))

row = 0

for per in all_per:

    if per == 'percent1':
        percent = 100
    elif per == 'percent0.95':
        percent = 95
    elif per == 'percent0.9':
        percent = 90
    elif per == 'percent0.8':
        percent = 80


    axs[row,0].set_xlabel(r'$x$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,0].set_xlim(minns[0], maxxs[0])

    axs[row,0].set_ylabel(r'$y$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,0].set_ylim(minns[1], maxxs[1])

    axs[row,1].set_xlabel(r'$z$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,1].set_xlim(minns[2], maxxs[2])

    axs[row,1].set_ylabel(r'$y$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,1].set_ylim(minns[1], maxxs[1])

    axs[row,2].set_xlabel(r'$x$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,2].set_xlim(minns[0], maxxs[0])

    axs[row,2].set_ylabel(r'$z$ (h$^{-1}$ Mpc)', fontsize=14)
    axs[row,2].set_ylim(minns[2], maxxs[2])

    
    for typ in sdss_info:

        this_typ = sdss_info[typ]

        this_per = this_typ[per]

        print(per, percent)

        idx2 = 0

        for sample in this_per:

            this_sample = this_per[sample]

            all_centers = []

            for boundary in this_sample['smoothened_minimal']:

                center = np.mean(boundary, axis=0)

                all_centers.append(center)

            all_centers = np.array(all_centers)

            axs[row,0].scatter(all_centers[:,0], all_centers[:,1]\
                    , alpha=0.5, marker=markers_dict[typ], color=colors_list[idx2])

            axs[row,1].scatter(all_centers[:,2], all_centers[:,1]\
                    , alpha=0.5, marker=markers_dict[typ], color=colors_list[idx2])

            axs[row,2].scatter(all_centers[:,0], all_centers[:,2]\
                    , alpha=0.5, marker=markers_dict[typ], color=colors_list[idx2])

            idx2 += 1

    custom_lines = [Line2D([0], [0], color='w', markeredgecolor='black'\
                            , markerfacecolor = 'none', marker='o', lw=4, alpha=0.8),
                    Line2D([0], [0], color='w', markeredgecolor='black'\
                            , marker='x', lw=4, alpha=0.8)]
    
    #axs[row,0].legend(custom_lines, ['Original', 'Trimmed'])

    row += 1

plt.tight_layout()

plt.savefig('figures/all_per_xyz_samethreshold.jpg', dpi=100)
plt.cla()
plt.clf()
plt.close()

# plot Mann-whitney
colors = dict()
markers = dict()

all_centers = dict()

markers_list = ['o', 'x']
colors_list = ['#D81B60', '#1E88E5', '#FFC107']

idx = 0
for typ in sdss_info:
    markers[typ] = markers_list[idx]
    idx2 = 0
    for sample in sdss_info[typ]['percent1']:
        colors[sample] = colors_list[idx2]
        idx2 += 1
    
    idx += 1

all_typ = list(sdss_info.keys())
all_per = list(sdss_info[all_typ[0]].keys())
all_sample = list(sdss_info[all_typ[0]][all_per[0]].keys())


all_cats = []

for per in all_per:
    #for per in ['percent0.9', 'percent0.8']:

    if per not in all_centers:
        all_centers[per] = dict()

    idx = 0

    for typ in all_typ:

        if typ not in all_centers[per]:
            all_centers[per][typ] = dict()


        this_typ = sdss_info[typ]

        this_per = this_typ[per]

        for sample in all_sample:

            this_sample = this_per[sample]

            all_cent = []
            
            for boundary in this_sample['smoothened_minimal']:

                center = np.mean(boundary, axis=0)
                
                all_cent.append(center)

            all_cent = np.array(all_cent)

            all_centers[per][typ][sample] = all_cent

        idx += 1

fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

axs[0].set_ylim([0, 1])
axs[1].set_ylim([0, 1])
axs[2].set_ylim([0, 1])

# Do original-original
typ = all_typ[0]
color = colors_list[0]

xx = 0
xx_tick = []
xx_label = []

per_labels = ['100', '95', '90', '85']

minn_p = math.inf

for idx1, p1 in enumerate(all_per):
    for idx2, p2 in enumerate(all_per):

        xx_tick.append(xx)

        
        xx_label.append(per_labels[idx1]+'-'+per_labels[idx2])

        all1 = all_centers[p1][typ]
        all2 = all_centers[p2][typ]


        all_this_typ_x = []
        all_this_typ_y = []
        all_this_typ_z = []

        if p1 != p2:
            ll = list(it.combinations_with_replacement(all_sample, 2))
        else:
            ll = list(it.combinations(all_sample, 2))

        for s1, s2 in ll:
            this_s1 = all1[s1]
            this_s2 = all2[s2]

            _, p = mannwhitneyu(this_s1, this_s2)

            if np.amin(p) < minn_p:
                minn_p = np.amin(p)

            all_this_typ_x.append(p[0])
            all_this_typ_y.append(p[1])
            all_this_typ_z.append(p[2])


        all_this_typ_x = np.array(all_this_typ_x)
        all_this_typ_y = np.array(all_this_typ_y)
        all_this_typ_z = np.array(all_this_typ_z)
        
        axs[0].errorbar(xx, y=np.mean(all_this_typ_x), yerr=np.std(all_this_typ_x)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[1].errorbar(xx, y=np.mean(all_this_typ_y), yerr=np.std(all_this_typ_y)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[2].errorbar(xx, y=np.mean(all_this_typ_z), yerr=np.std(all_this_typ_z)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        xx += 1


# Do trimmed-trimmed
typ = all_typ[1]
color = colors_list[1]
xx = 0
for p1 in all_per:
    for p2 in all_per:

        all1 = all_centers[p1][typ]
        all2 = all_centers[p2][typ]

        all_this_typ_x = []
        all_this_typ_y = []
        all_this_typ_z = []

        if p1 != p2:
            ll = list(it.combinations_with_replacement(all_sample, 2))
        else:
            ll = list(it.combinations(all_sample, 2))

        for s1, s2 in ll:
            this_s1 = all1[s1]
            this_s2 = all2[s2]

            _, p = mannwhitneyu(this_s1, this_s2)
            if np.amin(p) < minn_p:
                minn_p = np.amin(p)

            all_this_typ_x.append(p[0])
            all_this_typ_y.append(p[1])
            all_this_typ_z.append(p[2])


        all_this_typ_x = np.array(all_this_typ_x)
        all_this_typ_y = np.array(all_this_typ_y)
        all_this_typ_z = np.array(all_this_typ_z)
        
        axs[0].errorbar(xx, y=np.mean(all_this_typ_x), yerr=np.std(all_this_typ_x)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[1].errorbar(xx, y=np.mean(all_this_typ_y), yerr=np.std(all_this_typ_y)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[2].errorbar(xx, y=np.mean(all_this_typ_z), yerr=np.std(all_this_typ_z)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        xx += 1

# Do original-trimmed
xx = 0
color = colors_list[2]
for p1 in all_per:
    for p2 in all_per:

        all1 = all_centers[p1]
        all2 = all_centers[p2]

        all_this_typ_x = []
        all_this_typ_y = []
        all_this_typ_z = []

        for s1, s2 in list(it.combinations_with_replacement(all_sample, 2)):
            this_s1 = all1[all_typ[0]][s1]
            this_s2 = all2[all_typ[1]][s2]

            _, p = mannwhitneyu(this_s1, this_s2)
            if np.amin(p) < minn_p:
                minn_p = np.amin(p)

            all_this_typ_x.append(p[0])
            all_this_typ_y.append(p[1])
            all_this_typ_z.append(p[2])


        all_this_typ_x = np.array(all_this_typ_x)
        all_this_typ_y = np.array(all_this_typ_y)
        all_this_typ_z = np.array(all_this_typ_z)
        
        axs[0].errorbar(xx, y=np.mean(all_this_typ_x), yerr=np.std(all_this_typ_x)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[1].errorbar(xx, y=np.mean(all_this_typ_y), yerr=np.std(all_this_typ_y)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        axs[2].errorbar(xx, y=np.mean(all_this_typ_z), yerr=np.std(all_this_typ_z)\
                , fmt='o', color=color\
                , ecolor=color, elinewidth=3, capsize=0\
                , alpha=0.5)

        xx += 1



custom_lines = [Line2D([0], [0], color=colors_list[0], lw=4, marker = 'o', alpha=0.8),
                Line2D([0], [0], color=colors_list[1], lw=4, marker = 'o', alpha=0.8),
                Line2D([0], [0], color=colors_list[2], lw=4, marker = 'o', alpha=0.8)]

axs[0].legend(custom_lines, ['original-original', 'trimmed-trimmed', 'original-trimmed']\
        , loc='lower right')
axs[1].legend(custom_lines, ['original-original', 'trimmed-trimmed', 'original-trimmed']\
        , loc='lower right')
axs[2].legend(custom_lines, ['original-original', 'trimmed-trimmed', 'original-trimmed']\
        , loc='lower right')
axs[2].set_xticks(xx_tick, xx_label, rotation=30)
axs[2].set_xlabel('Pair of percentages compared', fontsize=14)
axs[0].set_ylabel('p-values', fontsize=14)
axs[1].set_ylabel('p-values', fontsize=14)
axs[2].set_ylabel('p-values', fontsize=14)

plt.savefig('figures/mannwhitney_all_samethreshold.jpg', dpi=600)






