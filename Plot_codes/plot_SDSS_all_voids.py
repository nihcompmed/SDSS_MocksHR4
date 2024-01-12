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
import helper_functions as hf


sdss_file = '../sdss_all_info_w_smoothened.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

tau_u = 22.5

# Sample1
locs = sdss_info['original_undersamples']['percent1']['sample1']['locs']
sdss_data = sdss_info['original_undersamples']['percent1']['sample1']['smoothened_minimal']

# Plot galaxies
mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None, engine=None, size=(400, 350))

norms = np.sum(locs**2,axis=-1)**(1./2)
#norms = norms + 1 

norm = mpl.colors.PowerNorm(gamma=1/3.75)
#norm = mpl.colors.Normalize()

cmap = cm.gnuplot
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
g.glyph.glyph.scale_factor = 1  # set scaling for all the points
g.glyph.scale_mode = 'data_scaling_off' # make all the points same size

## Visualize local density
#mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts))

count = 0

for this_boundary in sdss_data:

    index_dict = dict()

    pt_idx = 0

    bndry_pts = []

    for pt in this_boundary:
        if frozenset(pt) not in index_dict:
            index_dict[frozenset(pt)] = pt_idx
            bndry_pts.append(pt)
            pt_idx += 1

    print(len(this_boundary))
    assert len(this_boundary) % 3 == 0, 'what???'

    bndry_pts = np.array(bndry_pts)
    this_triangles = []
    nn = 0
    while nn < len(this_boundary):
        this_triangles.append([index_dict[frozenset(this_boundary[nn])]\
                               , index_dict[frozenset(this_boundary[nn+1])]\
                               , index_dict[frozenset(this_boundary[nn+2])]])
        nn += 3
    
    if count == 0:
        color = (75/255,170/255,210/255)
        opacity = 1
    else:
        color = (1,223/255,120/255)
        opacity = 0.65

    count += 1
    mlab.triangular_mesh(bndry_pts[:, 0]\
                        , bndry_pts[:, 1]\
                        , bndry_pts[:, 2]\
                        , this_triangles\
                        , color=color\
                        , opacity=opacity\
                        , representation='surface'\
                        )
#minimal_candidates.append(cover_locs[list(frozenset(cyc))])

mlab.points3d([0, 100], [0, 0], [0, 0], color=(1,0,0), scale_factor=10)

mlab.show()



