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



sdss_file = '../sdss_all_info_w_smoothened.p'
sdss_info = pickle.load(open(sdss_file, 'rb'))

tau_u = 22.5

# Sample1
locs = sdss_info['original_undersamples']['percent1']['sample1']['locs']
sdss_data = sdss_info['original_undersamples']['percent1']['sample1']['smoothened_minimal']


tree = scipy.spatial.KDTree(locs)

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

#print(f'Total grid points are {grid_points.shape[0]}')

# Plot galaxies
mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None, engine=None, size=(400, 350))

#norms = np.sum(locs**2,axis=-1)**(1./2)
#
#norm = mpl.colors.PowerNorm(gamma=1/3.75)
#
#cmap = cm.gnuplot
#m = cm.ScalarMappable(norm=norm, cmap=cmap)
#rgba = m.to_rgba(norms)
#
#rgba[:,-1] = 0.2
#
## Plot all points
#print('Plotting all points...')
#pts = mlab.pipeline.scalar_scatter(locs[:,0], locs[:,1], locs[:, 2]\
#                                , resolution=16) # plot the points
#
#pts.add_attribute(rgba, 'colors') # assign the colors to each point
#pts.data.point_data.set_active_scalars('colors')
#g = mlab.pipeline.glyph(pts, opacity=0.4)
#g.glyph.glyph.scale_factor = 2  # set scaling for all the points
#g.glyph.scale_mode = 'data_scaling_off' # make all the points same size


all_pts_inside = []

for this_boundary in sdss_data:

    pts_inside = grid_points[Delaunay(this_boundary).find_simplex(grid_points) >= 0]

    #all_pts_inside += list(pts_inside)

    #mlab.points3d(pts_inside[:,0]\
    #            , pts_inside[:,1]\
    #            , pts_inside[:,2]\
    #            , color=(1,1,1)\
    #            , scale_factor=0.002\
    #            , opacity=0.3\
    #            )

    ## nearest galaxy
    #dd, ii = tree.query(pts_inside, k=1)

    #print(dd)
    #exit()

    this_centroid = np.mean(pts_inside, axis=0)

    diff = (pts_inside - this_centroid)**2

    this_dist = np.sum((pts_inside - this_centroid)**2, axis=-1)**(1./2)


    #norms = np.sum(np.abs(locs)**2,axis=-1)**(1./2)
    #norms = norms + 1 
    
    norm = mpl.colors.PowerNorm(gamma=1/2.75)
    #norm = mpl.colors.Normalize()
    #norm = mpl.colors.LogNorm()
    
    #cmap = cm.gnuplot
    cmap2 = cm.hot_r
    #cmap2 = cmap
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap2)

    #
    rgba2 = m2.to_rgba(this_dist)

    #rgba2 = m2.to_rgba(dd)

    #rgba[:,-1] = 0.2
    
    # Plot all points
    print('Plotting all points...')
    pts2 = mlab.pipeline.scalar_scatter(pts_inside[:,0], pts_inside[:,1], pts_inside[:, 2]\
                                    , resolution=16) # plot the points

    pts2.add_attribute(rgba2, 'colors') # assign the colors to each point
    pts2.data.point_data.set_active_scalars('colors')
    g2 = mlab.pipeline.glyph(pts2, opacity=0.4)
    g2.glyph.glyph.scale_factor = 4 # set scaling for all the points
    g2.glyph.scale_mode = 'data_scaling_off' # make all the points same size

    #mlab.points3d([0, 100], [0, 0], [0, 0], scale_factor=10, color=(1,0,0))
    #mlab.show()




#ref_x = np.array([0,188])
#ref_y = np.array([0,0])
#ref_z = np.array([0,0])
#ref_pts = np.array([[0,0,0], [188,0,0]])
#tube = mlab.pipeline.tube(ref_pts, tube_radius=10)
#tube.filter.radius_factor = 1.
#tube.filter.vary_radius = 'vary_radius_by_scalar'
#mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

#mlab.pipeline.line_source(ref_x, ref_y, ref_z, color=(1,0,0))

mlab.points3d([0, 100], [0, 0], [0, 0], scale_factor=10, color=(1,0,0))
#mlab.points3d([188], [0], [0], scale_factor=0.01, color=(0,1,0))

#mlab.pipeline.tube(ref_x, ref_y, ref_z, tube_radius=1, tube_sides=12, color=(1,0,0))
mlab.show()



