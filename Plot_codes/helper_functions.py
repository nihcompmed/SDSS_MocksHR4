import numpy as np
from scipy.spatial import ConvexHull
import math
from scipy.spatial import Delaunay, KDTree
import networkx as nx
import itertools
from astropy.cosmology import FlatLambdaCDM
import astropy.cosmology.units as cu
import astropy.units as u
import astropy

def get_vol_convexhull(boundary):
    vol = ConvexHull(boundary).volume
    return vol

def get_reff_convexhull(boundary):
    vol = get_vol_convexhull(boundary)
    reff = (3*vol/(4*math.pi))**(1/3)
    return reff


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

# thresh is birth threshold
def smoothen_trivial_tetra(boundary, unique_pts, locs, thresh):

    # 1. Make a list of all triangles in boundary, format [{v1, v2, v3}]

    triangles = []

    idx = 0
    lenn = len(boundary)

    while idx < lenn:
        v1 = unique_pts[frozenset(boundary[idx])]
        v2 = unique_pts[frozenset(boundary[idx+1])]
        v3 = unique_pts[frozenset(boundary[idx+2])]

        this_triangle = frozenset({v1, v2, v3})

        #assert len(this_triangle) == 3, 'Triangle is not a triangle?'
        # Ignore degenerate triangles
        if len(this_triangle) == 3:
            # Ignore non-unique triangles
            if frozenset({v1, v2, v3}) not in triangles:
                triangles.append(frozenset({v1, v2, v3}))

        idx += 3


    # Make dual graph
    dual_graph = nx.Graph()

    # Add triangles as nodes
    for t in triangles:

        diam = 0
        for v1, v2 in itertools.combinations(list(t), 2):
            lenn = np.sum((locs[v1] - locs[v2])**2)
            if lenn > diam:
                diam = lenn

        diam = math.sqrt(diam)
        dual_graph.add_node(t, diam = diam)

    # Add adjacent triangles as edges of the dual graph
    for t1, t2 in itertools.combinations(triangles, 2):
        t3 = t1.intersection(t2)
        assert len(t3) < 3, 'how is it more than or equal to 3?'
        if len(t3) == 2:
            dual_graph.add_edge(t1, t2)


    update = 1
    while (update):

        #####################
        #### FOR TESTING

        ###plot_triangles = []
        ###for node in dual_graph.nodes:
        ###    plot_triangles.append(list(node))

        ###mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
        ###                        , engine=None, size=(400, 350))
        ###
        ###mlab.points3d(locs[:,0]\
        ###                        , locs[:,1]\
        ###                        , locs[:,2]\
        ###                        , scale_factor = 0.2)

        ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
        ###                , plot_triangles, color=(1,223/255,120/255), opacity=0.3)

        ###mlab.show()

        #####################


        update = 0

        # Get all tetraherons, a clique of t1, t2, t3, t4

        # Fist, get all maximal cliques
        #print('Finding cliques')
        cliques = nx.find_cliques(dual_graph)

        # Remove the cliques that have less than 3 nodes
        useful_cliques = [x for x in cliques if len(x) > 2]

        smallest_diam = math.inf
        small_tetra = []
        
        ## IGNORE TETRA REMOVAL
        ## Now, go over cliques of len > 3 to get tetrahedrons
        #for cl in useful_cliques:
        #    if len(cl) < 4:
        #        continue

        #    for t1, t2, t3, t4 in itertools.combinations(list(cl), 4):
        #        diam = min(dual_graph.nodes[t1]['diam']\
        #                , dual_graph.nodes[t2]['diam']\
        #                , dual_graph.nodes[t3]['diam']\
        #                , dual_graph.nodes[t4]['diam']\
        #                )
        #        if diam < smallest_diam:
        #            smallest_diam = diam
        #            small_tetra = [t1, t2, t3,t4]

        ## if there was a smallest tetra, remove it from dual graph
        #if smallest_diam != math.inf:
        #    dual_graph.remove_nodes_from(small_tetra)
        #    update = 1
        #    continue

        # No tetra found in boundary, so, look for triangles now
        for cl in useful_cliques:
            for t1, t2, t3 in itertools.combinations(list(cl), 3):

                # This should not be a tetra
                # So, check if new_face is already in the clique or not
                all_verts = t1.union(t2)
                all_verts = all_verts.union(t3)

                new_face = all_verts - t1
                new_face = new_face.union(all_verts - t2)
                new_face = new_face.union(all_verts - t3)

                if new_face in list(cl):
                    # This means tetrahedron [t1, t2, t3, new_face] is in the boundary. Skip it.
                    continue

                assert len(new_face) == 3, 'new face is not a triangle???'
                
                # new_face is not in the boundary

                # Find diameter of new_face:
                new_face_diam = 0
                for v1, v2 in itertools.combinations(list(new_face), 2):
                    lenn = np.sum((locs[v1] - locs[v2])**2)
                    if lenn > new_face_diam:
                        new_face_diam = lenn
                new_face_diam = math.sqrt(new_face_diam)

                # Diameter of the tetrahedron is the largest of diameter of its faces
                # (Alternately, could have looked at lengths of all edges in the tetrahedron)
                diam = max(dual_graph.nodes[t1]['diam']\
                        , dual_graph.nodes[t2]['diam']\
                        , dual_graph.nodes[t3]['diam']\
                        , new_face_diam\
                        )

                # Consider only if the diameter of this tetrahedron is <= thresh
                if (diam < smallest_diam) and (diam <= thresh):
                    smallest_diam = diam
                    small_tetra = [t1, t2, t3]
                    small_new_face = new_face
                    small_new_face_diam = new_face_diam

        # If there was a smallest triangle with diameter, then remove t1,t2,t3 from dual, BUT...
        # ...have to add the fourth face/triangle
        if smallest_diam != math.inf:

            t1, t2, t3 = small_tetra

            new_face = small_new_face

            # Add the new face
            #print('adding', new_face)
            dual_graph.add_node(new_face, diam = small_new_face_diam)

            #mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
            #                        , engine=None, size=(400, 350))

            #mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            #        ,[list(new_face)], color=(0,0,1), opacity=0.8)

            # Go over neighbors of t1, t2, t3
            neigh = list(dual_graph.neighbors(t1))
            neigh += list(dual_graph.neighbors(t2))
            neigh += list(dual_graph.neighbors(t3))
            neigh = frozenset(neigh)

            #print('t1', t1)
            #print('t2', t2)
            #print('t3', t3)
            #print('new_face', new_face)

            for n in neigh:
                kk = n.intersection(new_face)
                #print(new_face, n, kk)
                #assert len(kk) < 3, 'what??????'
                if len(kk) == 2:
                    dual_graph.add_edge(new_face, n)

            # Remove t1, t2, t3
            #print('removing', t1, t2, t3)
            dual_graph.remove_nodes_from([t1, t2, t3])

            
            #######################
            ##### For testing
            ###triangles = []
            ###triangles.append(list(t1))
            ###triangles.append(list(t2))
            ###triangles.append(list(t3))

            ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            ###        ,triangles, color=(1,0,0), opacity=0.3)

            ###
            ###mlab.points3d(locs[:,0]\
            ###                        , locs[:,1]\
            ###                        , locs[:,2]\
            ###                        , scale_factor = 0.2)

            ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            ###                , plot_triangles, color=(1,223/255,120/255), opacity=0.3)
            ###mlab.show()

            update = 1
            continue


    # Remove 'disconnected' tetrahedrons

    dual_new = []

    for comp in nx.connected_components(dual_graph):
        # If there are <= 4 triangle -> there is no non-trivial boundary
        if len(comp) <= 4:
            #print('skipping these points')
            continue
        dual_new += comp

    dual_graph = dual_graph.subgraph(dual_new)


    new_boundary_list = []

    for node in dual_graph.nodes:
        new_boundary_list += [locs[n] for n in list(node)]

    
    #mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
    #                        , engine=None, size=(400, 350))
    
    #mlab.points3d(locs[:,0]\
    #                        , locs[:,1]\
    #                        , locs[:,2]\
    #                        , scale_factor = 0.2)
    
    #plot_H2_boundary(locs, new_boundary_list, mlab, (238/255, 75/255, 43/255), 0.5)


    #mlab.show()

    #exit()

    return np.array(new_boundary_list)


def local_smoothen_boundary(boundary, thresh):

        index_pts = dict()
        unique_pts = dict()

        idx = 0

        for pt in boundary:
            
            set_pt = frozenset(pt)
            if set_pt not in unique_pts:
                unique_pts[set_pt] = idx
                index_pts[idx] = pt
                idx += 1

        new_boundary = smoothen_trivial_tetra(boundary, unique_pts, index_pts, thresh)

        return new_boundary



# Cartesian to ra, dec, redshift
def invert(xx, yy, zz, cosmo):

    #cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    val = xx**2 + yy**2

    z = math.sqrt(val + zz**2)

    d = math.asin(zz/z)
    d = d*180/math.pi

    a = abs(math.atan(yy/xx))
    a = a*180/math.pi

    if xx < 0 and yy > 0:
        a = 180 - a
    elif xx < 0 and yy < 0:
        a = 180 + a
    elif xx > 0 and yy < 0:
        a = 360 - a

    z = z * u.Mpc

    z = z.to(cu.redshift, cu.redshift_distance(cosmo, kind="comoving"))

    return a, d, z.value



