from mayavi import mlab
import pickle
import numpy as np
import matplotlib.pyplot as plt
import umap
from scipy.stats import mannwhitneyu
import itertools as it
import math
import seaborn as sns
from matplotlib.lines import Line2D
from astropy.cosmology import FlatLambdaCDM
import astropy.cosmology.units as cu
import astropy.units as u
import astropy
from scipy.spatial import ConvexHull


#print(astropy.__version__)
cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
#print(cosmo)



def invert(xx, yy, zz):

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

    return a, d, z



plt.rcParams['text.usetex'] = True


all_info = pickle.load(open('../samethreshold_all_info_w_smoothened.p', 'rb'))

#reducer = umap.UMAP(random_state=42)


all_centers = dict()


all_typ = list(all_info.keys())
all_per = list(all_info[all_typ[0]].keys())
all_sample = list(all_info[all_typ[0]][all_per[0]].keys())

#fig, axs = plt.subplots(4, 2, figsize=(10,8), sharex=True, sharey=True)



for per in all_per:

    if per != 'percent1':
        continue


    for typ in all_info:
        print(typ)

        this_typ = all_info[typ]

        this_per = this_typ[per]

        idx2 = 0

        for sample in this_per:

            this_sample = this_per[sample]

            all_centers = []

            for boundary in this_sample['smoothened_minimal']:

                vol = ConvexHull(boundary).volume

                eff_rad = (3*vol/(4*math.pi))**(1/3)

                center = np.mean(boundary, axis=0)

                a, d, z = invert(center[0], center[1], center[2])

                all_centers.append([a, d, z , eff_rad])

            all_centers = np.array(all_centers)

            all_centers = all_centers[all_centers[:, 0].argsort()]

            fname = typ + '_' + sample + '.csv' 

            np.savetxt(fname, all_centers, delimiter='&', fmt='%.4f')


    

                    
