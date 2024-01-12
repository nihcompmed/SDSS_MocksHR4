import numpy as np
import os
from mayavi import mlab
import matplotlib.cm as cm
import matplotlib as mpl
import pickle
from scipy.spatial import Delaunay, KDTree
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import scipy
import seaborn as sns
from sklearn.metrics import pairwise_distances
from joblib import Parallel, delayed
#import smoothen_tetra as smoothen
from scipy.spatial.distance import cdist as sci_cdist
from scipy.spatial import ConvexHull
import math
import itertools
import networkx as nx
from tqdm import tqdm
import helper_functions as hf



percent_vol_overlap_thresh = 2

sdss_file = '../sdss_all_info_w_smoothened.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))


# Sample1
locs = sdss_info['original_undersamples']['percent1']['sample1']['locs']
sdss_data = sdss_info['original_undersamples']['percent1']['sample1']['smoothened_minimal']

# hyper rectangle
minns = np.min(locs, axis = 0)
maxxs = np.max(locs, axis = 0)

grid_reso = 3

n_points = np.int32((maxxs - minns)/grid_reso) + 1

xx = np.linspace(minns[0], maxxs[0], n_points[0])
yy = np.linspace(minns[1], maxxs[1], n_points[1])
zz = np.linspace(minns[2], maxxs[2], n_points[2])

xv, yv, zv = np.meshgrid(xx, yy, zz)

grid_points = np.vstack((xv.flatten(), yv.flatten(), zv.flatten())).T

grid_tree = KDTree(grid_points)


# Plot galaxies
mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None, engine=None, size=(400, 350))

norms = np.sum(np.abs(locs)**2,axis=-1)**(1./2)
#norms = norms + 1 

norm = mpl.colors.PowerNorm(gamma=1/3.75)
#norm = mpl.colors.Normalize()

cmap = cm.gnuplot
##cmap = cm.Greys
m = cm.ScalarMappable(norm=norm, cmap=cmap)
#m = cm.ScalarMappable(cmap=cmap)
rgba = m.to_rgba(norms)

rgba[:,-1] = 0.2

# Plot all points
print('Plotting all points...')
pts = mlab.pipeline.scalar_scatter(locs[:,0], locs[:,1], locs[:, 2]\
                                , resolution=16) # plot the points

pts.add_attribute(rgba, 'colors') # assign the colors to each point
pts.data.point_data.set_active_scalars('colors')
g = mlab.pipeline.glyph(pts, opacity=0.4)
g.glyph.glyph.scale_factor = 2 # set scaling for all the points
g.glyph.scale_mode = 'data_scaling_off' # make all the points same size

PH_centers = []
PH_radii = []
PH_boundary = []

for this_boundary in tqdm(sdss_data):

    this_center = np.mean(this_boundary, axis=0)
    
    PH_centers.append(this_center)
    PH_boundary.append(this_boundary)


PH_centers = np.array(PH_centers)



V2_info_file = '../V2_analysis/V2_VIDE/processed_info.p'
V2_info = pickle.load(open(V2_info_file, 'rb'))


V2_dirr = '../V2_analysis/V2_voids'
n_V2 = 419

shortlisted_V2 = []

# Thresholds in h^-1 Mpc
epsilon = 7.5
tau_u = 22.5

count = 0
skip_count = 0
V2_sig_count = 0

for v_idx in range(n_V2):
    void_dirr = V2_dirr + '/void'+str(v_idx)
    if not os.path.isdir(void_dirr):
        continue

    locs_file = void_dirr + '/locs.csv'
    if not os.path.isfile(locs_file):
        skip_count += 1

    pers_file = void_dirr + '/H2_pers_data.txt'
    if not os.path.isfile(pers_file):
        continue

    pers = np.loadtxt(pers_file, delimiter=',')
    if not len(pers):
        continue


    if pers.ndim == 1:
        pers = np.reshape(pers, (1, 2))

    pers[pers[:,1]==-1, 1] = epsilon + tau_u

    diff = pers[:,1] - pers[:,0]

    sig_idxs = np.argwhere(diff >= epsilon).flatten()

    n_sig = len(sig_idxs)

    if not n_sig:
        continue

    V2_sig_count += 1

    boundary_pts = []
    for triangle in V2_info['boundary'][v_idx]:
        for pt in triangle:
            #if pt not in boundary_pts:
            boundary_pts.append(pt)

    # THERE HAS TO BE A FASTER WAY TO DO THIS USING NUMPY

    shortlisted_V2.append(boundary_pts)
    count += 1
    

print(f'{skip_count} voids for which V2 did not report a boundary')
print(f'{V2_sig_count} V2 voids  contain at least one sig feature')

sig_V2_centers = []
for cyc in shortlisted_V2:
    this_center = np.mean(cyc, axis=0)
    sig_V2_centers.append(this_center)

sig_V2_centers = np.array(sig_V2_centers)

# Volume inside PH
all_PH_pts_inside = []
all_PH_idxs_inside = []
all_PH_galaxies_inside = []
all_PH_vol = []

all_PH_gor = []

for boundary in PH_boundary:

    gor_r, PH_gor = hf.get_gor(boundary, locs)

    all_PH_gor.append(PH_gor)


    this_idxs_inside = np.argwhere(Delaunay(boundary).find_simplex(grid_points) >= 0).flatten()
    all_PH_idxs_inside.append(list(this_idxs_inside))

    this_gal_inside = np.argwhere(Delaunay(boundary).find_simplex(locs) >= 0).flatten()
    all_PH_galaxies_inside.append(len(this_gal_inside))

    pts_inside = grid_points[Delaunay(boundary).find_simplex(grid_points) >= 0]
    
    all_PH_pts_inside += list(pts_inside)

    all_PH_vol.append(hf.get_vol_convexhull(boundary))

all_PH_pts_inside = np.array(all_PH_pts_inside)


mlab.points3d(all_PH_pts_inside[:,0]\
            , all_PH_pts_inside[:,1]\
            , all_PH_pts_inside[:,2]\
            , opacity=0.05\
            , scale_factor=5\
            )

# Volume inside V2
all_V2_pts_inside = []
all_V2_idxs_inside = []
all_V2_galaxies_inside = []
all_V2_vol = []

all_V2_gor = []

for boundary in shortlisted_V2:

    gor_r, V2_gor = hf.get_gor(boundary, locs)

    all_V2_gor.append(V2_gor)

    this_idxs_inside = np.argwhere(Delaunay(boundary).find_simplex(grid_points) >= 0).flatten()
    all_V2_idxs_inside.append(list(this_idxs_inside))

    this_gal_inside = np.argwhere(Delaunay(boundary).find_simplex(locs) >= 0).flatten()
    all_V2_galaxies_inside.append(len(this_gal_inside))

    pts_inside = grid_points[Delaunay(boundary).find_simplex(grid_points) >= 0]
    
    all_V2_pts_inside += list(pts_inside)

    all_V2_vol.append(hf.get_vol_convexhull(boundary))

all_V2_pts_inside = np.array(all_V2_pts_inside)

mlab.points3d(all_V2_pts_inside[:,0]\
            , all_V2_pts_inside[:,1]\
            , all_V2_pts_inside[:,2]\
            , opacity=0.05\
            , scale_factor=5\
            , color=(0,0,1)\
            )

mlab.points3d([0, 100], [0, 0], [0, 0], scale_factor=10, color=(1,0,0))
#mlab.show()
#exit()

all_PH_gor = np.array(all_PH_gor)
all_V2_gor = np.array(all_V2_gor)

for this_gor in all_PH_gor:
    plt.plot(gor_r, this_gor, alpha=0.2, color='red')

gor_PH_med = np.median(all_PH_gor, axis=0)
plt.plot(gor_r, gor_PH_med, color='red', lw=2, label='PH voids median')

for this_gor in all_V2_gor:
    plt.plot(gor_r, this_gor, alpha=0.2, color='blue')

gor_V2_med = np.median(all_V2_gor, axis=0)
plt.plot(gor_r, gor_V2_med, color='blue', lw=2, label='V$^2$ voids median')

plt.xlabel(r'$r$ ($h^{-1}$ Mpc)', fontsize=16)
plt.ylabel(r'$g(r)$', fontsize=16)

plt.legend()

plt.savefig(f'figures/SDSS_sample1_V2_gor.jpg', dpi=600)
plt.cla()
plt.clf()

# Mean with std. dev

PH_mean = np.mean(all_PH_gor, axis=0)
PH_std = np.std(all_PH_gor, axis=0)

plt.plot(gor_r, PH_mean, color='red', lw=2, label=r'PH voids mean')
plt.fill_between(gor_r, PH_mean-PH_std, PH_mean+PH_std\
                , zorder=2000, alpha=0.4, color='red', label=r'PH std. dev.')

V2_mean = np.mean(all_V2_gor, axis=0)
V2_std = np.std(all_V2_gor, axis=0)
plt.plot(gor_r, V2_mean, color='blue', lw=2, label=r'V$^2$ voids mean')
plt.fill_between(gor_r, V2_mean-V2_std, V2_mean+V2_std\
                , zorder=2000, alpha=0.4, color='blue', label=r'V$^2$ std. dev.')

plt.xlabel(r'$r$ ($h^{-1}$ Mpc)', fontsize=18)
plt.ylabel(r'$g(r)$', fontsize=18)
plt.tight_layout()
plt.legend()
plt.savefig(f'figures/SDSS_sample1_V2_gor_mean_std.jpg', dpi=600)
plt.cla()
plt.clf()

#PH_pts_set = frozenset(all_PH_idxs_inside)
#V2_pts_set = frozenset(all_V2_idxs_inside)
#common = PH_pts_set.intersection(V2_pts_set)
#
#print(len(common)/len(PH_pts_set), len(common)/len(V2_pts_set))

#print(count)

PH_info_dict = dict()
V2_info_dict = dict()

PH_info_vol_match_dict = dict()
V2_info_vol_match_dict = dict()

G = nx.Graph()

for (b1,b2) in itertools.product(all_PH_idxs_inside, all_V2_idxs_inside):
    b1 = frozenset(b1)
    b2 = frozenset(b2)

    common = b1.intersection(b2)

    b1_vol = len(b1)
    b2_vol = len(b2)
    co_vol = len(common)

    b1_percent = co_vol/b1_vol*100
    b2_percent = co_vol/b2_vol*100

    min_percent = min(b1_percent, b2_percent)
    if min_percent >= 75:
        G.add_edge(b1, b2)



    if b1 not in PH_info_dict:
        PH_info_dict[b1] = b1_percent
        PH_info_vol_match_dict[b1] = b2_percent
    if b2 not in V2_info_dict:
        V2_info_dict[b2] = b2_percent
        V2_info_vol_match_dict[b2] = b1_percent

    if b1_percent > PH_info_dict[b1]:
        PH_info_dict[b1] = b1_percent
        PH_info_vol_match_dict[b1] = b2_percent
    if b2_percent > V2_info_dict[b2]:
        V2_info_dict[b2] = b2_percent
        V2_info_vol_match_dict[b2] = b1_percent

    #print(b1_percent, b2_percent)

#nx.draw(G)
#plt.show()
#plt.cla()
#plt.clf()

print('two way vol intersectin percent for PH')
xx = []
yy = []

xx_special = []
yy_special = []

for key in PH_info_dict:
    this_x = PH_info_dict[key]
    this_y = PH_info_vol_match_dict[key]
    if this_x < percent_vol_overlap_thresh and this_y < percent_vol_overlap_thresh:
        #print(xx[-1], yy[-1])
        xx_special.append(this_x)
        yy_special.append(this_y)
    else:
        xx.append(this_x)
        yy.append(this_y)

plt.scatter(xx, yy, color='red', label = 'PH voids', alpha=0.5)
plt.plot(xx_special, yy_special, color='red', marker='*', alpha=0.5, markersize=14, lw=0)

print('two way vol intersection percent for V2')
xx = []
yy = []
for key in V2_info_dict:
    xx.append(V2_info_dict[key])
    yy.append(V2_info_vol_match_dict[key])
    if xx[-1] < percent_vol_overlap_thresh and yy[-1] < percent_vol_overlap_thresh:
        print(xx[-1], yy[-1])

plt.scatter(xx, yy, color='blue', label = r'V$^2$ voids', alpha=0.5)

plt.xlabel('Maximal % volume intersection wrt void', fontsize=14)
plt.ylabel('Maximal % volume intersection wrt the intersecting void', fontsize=10)

plt.legend()

#plt.show()
plt.savefig('figures/max_vol_intersections_twoway.pdf', dpi=600)
plt.cla()
plt.clf()
#for key in PH_info_dict:
#    print(PH_info_dict[key])
#
#for key in V2_info_dict:
#    print(V2_info_dict[key])

PH_vals = np.array(list(PH_info_dict.values()))
V2_vals = np.array(list(V2_info_dict.values()))

print('min PH overlap', np.amin(PH_vals))
print('min V2 overlap', np.amin(V2_vals))

print('max PH overlap', np.amax(PH_vals))
print('max V2 overlap', np.amax(V2_vals))


small_PH_overlap = PH_vals[PH_vals<percent_vol_overlap_thresh]
small_V2_overlap = V2_vals[V2_vals<percent_vol_overlap_thresh]


print(f'PH with less than {percent_vol_overlap_thresh}% overlap', len(small_PH_overlap))
print(f'V2 with less than {percent_vol_overlap_thresh}% overlap', len(small_V2_overlap))

large_PH_overlap = PH_vals[PH_vals>80]
large_V2_overlap = V2_vals[V2_vals>80]

print('PH with more than 80% overlap', len(large_PH_overlap))
print('V2 with more than 80% overlap', len(large_V2_overlap))

plt.scatter([0]*len(PH_info_dict), list(PH_info_dict.values()), color='red', alpha=0.5)
plt.scatter([1]*len(V2_info_dict), list(V2_info_dict.values()), color='blue', alpha=0.5)
plt.ylabel('Maximum percentage volume intersection with other kind')

plt.xticks(ticks=[0, 1], labels=['PH', r'V$^2$'], fontsize=16)

plt.xlim([-0.5,1.5])

#plt.savefig('figures/vol_intersect_PH_V2.pdf', dpi=600)
plt.cla()
plt.clf()


#plt.scatter(all_PH_galaxies_inside, np.array(all_PH_galaxies_inside)/np.array(all_PH_vol), color='red', alpha=0.5, label='PH void')
#plt.scatter(all_V2_galaxies_inside, np.array(all_V2_galaxies_inside)/np.array(all_V2_vol), color='blue', alpha=0.5, label='V2 void')
#
#plt.show()

# galaxy containment box plot
cats = []
vals = []

PH_dens = np.array(all_PH_galaxies_inside)/np.array(all_PH_vol)
vals += list(PH_dens)
cats += ['PH']*len(all_PH_galaxies_inside)

V2_dens = np.array(all_V2_galaxies_inside)/np.array(all_V2_vol)
vals += list(V2_dens)
cats += [r'V$^2$']*len(all_V2_galaxies_inside)

print(scipy.stats.mannwhitneyu(PH_dens, V2_dens))

colors = dict()
colors['PH'] = 'red'
colors[r'V$^2$'] = 'turquoise'
sns.boxplot(x=cats, y=vals, palette=colors)
plt.gca().tick_params(axis='x', labelsize=16)

sns.stripplot(x=cats, y=vals, color='black')
plt.ylabel('Density (galaxies inside/volume)', fontsize=16)
#plt.savefig('figures/gal_inside_PH_V2.pdf', dpi=600)
plt.savefig('figures/gal_density_PH_V2.pdf', dpi=600)
plt.cla()
plt.clf()



#plt.show()
#exit()

#plt.scatter([0]*len(all_PH_galaxies_inside), all_PH_galaxies_inside, color='black', alpha=0.5)
#plt.scatter([1]*len(all_V2_galaxies_inside), all_V2_galaxies_inside, color='blue', alpha=0.5)

#plt.xticks(ticks=[0, 1], labels=['PH', r'V$^2$'], fontsize=16)

#plt.xlim([-0.5,1.5])

#plt.savefig('galaxies_inside_PH_V2.pdf', dpi=600)
#plt.cla()
#plt.clf()



#plt.show()






# Load V2 voids



#PH_radii = np.array(PH_radii)

#matched_indices = []

#print(V2_centers)

#mlab.points3d(V2_centers[:,0]\
#            , V2_centers[:,1]\
#            , V2_centers[:,2]\
#            , scale_factor=10\
#            , color=(67/255,232/255,216/255)\
#            , opacity=0.5\
#            )


#ddist = sci_cdist(PH_centers, V2_centers)

#print(ddist.shape)

### TRY MATCHING
#true_center_nearest = np.argmin(ddist, axis=0)
#estimated_center_nearest = np.argmin(ddist, axis=1)
#
#count = 0
#
#matched_V2_rad = []
#matched_PH_rad = []
#matched_dist = []
#
#for idx1, idx2 in enumerate(true_center_nearest):
#    if estimated_center_nearest[idx2] == idx1:
#        mlab.plot3d([V2_centers[idx1, 0], PH_centers[idx2, 0]]\
#                , [V2_centers[idx1, 1], PH_centers[idx2, 1]]\
#                , [V2_centers[idx1, 2], PH_centers[idx2, 2]]\
#                , color=(64/255,224/255,208/255), tube_radius=1)
#
#        count += 1
#        matched_V2_rad.append(V2_radii[idx1])
#        matched_PH_rad.append(PH_radii[idx2])
#
#        matched_indices.append(idx1)
#
#        matched_dist.append(math.sqrt(np.sum((V2_centers[idx1] - PH_centers[idx2])**2)))
#
#        #print(V2_radii[idx1], PH_radii[idx2], matched_dist)

#print('matched', count)

#print(matched_dist)

#unmatched_nsig = []
#matched_nsig = []
#
#count_n_sig = 0
#
#sig_centers = []
#insig_centers = []
#
#sig_idx = []
#
#for idx in range(152):
#
#
#    dirr = 'V2_void'+str(idx)
#
#    this_center = V2_centers[idx]
#
#    pers_file = dirr + '/H2_pers_data.txt'
#
#    pers = np.loadtxt(pers_file, delimiter=',')
#
#    pers[pers[:,1]==-1, 1] = epsilon + tau_u
#
#    diff = pers[:,1] - pers[:,0]
#
#    sig = np.argwhere(diff >= epsilon).flatten()
#
#    n_sig = len(sig)
#
#    if n_sig:
#        count_n_sig += n_sig
#        sig_centers.append(this_center)
#        sig_idx.append(idx)
#    else:
#        insig_centers.append(this_center)
#
#    #if idx in matched_indices:
#    #    matched_nsig.append(n_sig)
#    #else:
#    #    unmatched_nsig.append(n_sig)
#
##print(matched_nsig)
##print(unmatched_nsig)
##print(count_n_sig)
#
#sig_centers = np.array(sig_centers)
#insig_centers = np.array(insig_centers)
#
#mlab.points3d(sig_centers[:,0]\
#            , sig_centers[:,1]\
#            , sig_centers[:,2]\
#            , scale_factor=10\
#            , color=(255/255,0/255,0/255)\
#            , opacity=0.8\
#            )
#
#mlab.points3d(insig_centers[:,0]\
#            , insig_centers[:,1]\
#            , insig_centers[:,2]\
#            , scale_factor=10\
#            , color=(67/255,232/255,216/255)\
#            , opacity=0.4\
#            )
#
#all_PH_pts_inside = []
#
#for idx in range(len(PH_centers)):
#
#    this_PH_center = PH_centers[idx]
#
#    this_PH_reff = PH_radii[idx]
#
#    pts_inside = grid_tree.query_ball_point(this_PH_center, r=this_PH_reff)
#
#    all_PH_pts_inside += pts_inside
#
#    all_PH_pts_inside = list(frozenset(all_PH_pts_inside))
#
#all_PH_pts_inside_locs = grid_points[all_PH_pts_inside]
#
#mlab.points3d(all_PH_pts_inside_locs[:,0]\
#            , all_PH_pts_inside_locs[:,1]\
#            , all_PH_pts_inside_locs[:,2]\
#            , opacity=0.1\
#            , scale_factor=5\
#            )
#
#
#all_V2_pts_inside = []
#
#for idx in range(len(V2_centers)):
#
#    if idx not in sig_idx:
#        continue
#
#    this_V2_center = V2_centers[idx]
#
#    this_V2_max_rad = V2_max_radii[idx]
#    this_V2_reff = V2_reff[idx]
#
#    pts_inside = grid_tree.query_ball_point(this_V2_center, r=this_V2_reff)
#    #pts_inside = grid_tree.query_ball_point(this_V2_center, r=this_V2_max_rad)
#
#    all_V2_pts_inside += pts_inside
#
#    all_V2_pts_inside = list(frozenset(all_V2_pts_inside))
#
#all_V2_pts_inside_locs = grid_points[all_V2_pts_inside]
#
##PH_pts_set = frozenset(all_PH_pts_inside)
##V2_pts_set = frozenset(all_V2_pts_inside)
##
##common = PH_pts_set.intersection(V2_pts_set)
#
##print(len(common)/len(all_PH_pts_inside), len(common)/len(all_V2_pts_inside))
##exit()
#
#mlab.points3d(all_V2_pts_inside_locs[:,0]\
#            , all_V2_pts_inside_locs[:,1]\
#            , all_V2_pts_inside_locs[:,2]\
#            , opacity=0.1\
#            , scale_factor=5\
#            , color=(0,0,1)\
#            )
#
#mlab.show()
#
#exit()
#
#
#
#matched_V2_rad = np.array(matched_V2_rad)
#matched_PH_rad = np.array(matched_PH_rad)
#
#minn = min(np.amin(matched_V2_rad), np.amin(matched_PH_rad))
#maxx = max(np.amax(matched_V2_rad), np.amax(matched_PH_rad))
#
#
##plt.scatter(matched_V2_rad, matched_PH_rad)
#
#
##plt.show()
#sns.scatterplot(x=matched_V2_rad, y=matched_PH_rad, size=matched_dist)
#
#plt.plot([minn, maxx], [minn, maxx], ls='--')
#
#plt.xlabel('V2/VIDE Reff', fontsize=14)
#plt.ylabel('PH Reff', fontsize=14)
#
#plt.savefig('compare_V2_PH_Reff.pdf', dpi=600)
#
##plt.show()
#
#
