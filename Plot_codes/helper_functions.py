import numpy as np
from scipy.spatial import ConvexHull
import math
from scipy.spatial import Delaunay, KDTree



def get_reff_convexhull(boundary):

    vol = ConvexHull(boundary).volume
    return (3*vol/(4*math.pi))**(1/3)


def get_delaunay(points, points_to_check):
    Delaunay(tri).find_simplex(point)

def get_reff_delaunay(boundary):
    minn = np.min(boundary, axis=0)
    maxx = np.max(boundary, axis=0)

    x_ = np.arange(minn[0], maxx[0], 1)
    y_ = np.arange(minn[1], maxx[1], 1)
    z_ = np.arange(minn[2], maxx[2], 1)

    pts = []

    for x in x_:
        for y in y_:
            for z in z_:
                pts.append([x, y, z])

    pts = np.array(pts)

    test = Delaunay(boundary).find_simplex(pts) >= 0
    vol = test.sum()
    return (3*vol/(4*math.pi))**(1/3)

def get_gor(boundary, locs):

    locs_tree = KDTree(locs)

    center = np.mean(boundary, axis=0)

    # Get distance from center to farthest point
    #dists = scipy.spatial.distance.cdist(np.reshape(center, (1,3)), sdss_locs)
    #dists = scipy.spatial.distance.cdist(np.reshape(center, (1,3)), cyc_locs)
    #max_dist = np.amax(dists)+50
    #min_dist = np.amin(dists)

    RR = 100
    initial_R = 0

    dr = 2

    res = locs_tree.query_ball_point(np.reshape(center, (1,3)), r=initial_R)[0]
    count_prev = len(res)

    vol_prev = 4/3*math.pi*(initial_R**3)

    all_dens = []

    r_vals = np.arange(initial_R+dr, RR + dr, step=dr)

    N_all = len(locs_tree.query_ball_point(np.reshape(center, (1,3)), r=RR+dr)[0])

    rho = N_all/((4/3*math.pi*((RR+dr)**3)))

    for rr in r_vals:

        #print(rr, end='\r')

        # Neighbors in sphere of radius rr
        res = locs_tree.query_ball_point(np.reshape(center, (1,3)), r=rr)[0]
        count_next = len(res)
        vol_next = 4/3*math.pi*(rr**3)

        n_this = count_next - count_prev

        this_dens = n_this/(vol_next - vol_prev)

        this_dens = this_dens/rho

        all_dens.append(this_dens)

        count_prev = count_next
        vol_prev = vol_next


    return r_vals, np.array(all_dens)


def get_nearest_void_center_distances(boundaries):

    centroids = []

    for boundary in boundaries:
        centroids.append(np.mean(boundary, axis=0))

    centroids = np.array(centroids)

    centroids_tree = KDTree(centroids)

    dists, ind = centroids_tree.query(centroids, k=2)

    dists = dists[:,1]

    return dists







