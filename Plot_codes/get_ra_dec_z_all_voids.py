import numpy as np
import os
import math
import pickle
import helper_functions as hf
from tqdm import tqdm
from astropy.cosmology import FlatLambdaCDM
import astropy.cosmology.units as cu
import astropy.units as u
import astropy

sdss_file = '../sdss_all_info_w_smoothened.p'
all_info = pickle.load(open(sdss_file, 'rb'))

H0 = 100
Om0 = 0.3
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

gg = open('all_centroid_reff.csv', 'w')

for sample in ['sample0', 'sample1', 'sample2']:

    sdss_data = all_info['original_undersamples']['percent1'][sample]['smoothened_minimal']

    void_idx = 1

    gg.write(f'#SDSS_{sample}\n')
    gg.write(f'#void_idx,ra,dec,z,Reff\n')


    for this_boundary in sdss_data:

        #ff = open(f'SDSS_{sample}_void{void_idx}.csv', 'w')
        #ff.write(f'#ra,dec,z\n')

        #for pt in this_boundary:
        #    ra, dec, z = hf.invert(pt[0], pt[1], pt[2], cosmo)
        #    ll = ','.join([str(ra), str(dec), str(z)])
        #    ff.write(ll+'\n')

        #ff.close()

        centroid = np.mean(this_boundary, axis=0)
        
        xx = centroid[0]
        yy = centroid[1]
        zz = centroid[2]

        this_reff = hf.get_reff_convexhull(this_boundary)

        ra, dec, z = hf.invert(xx, yy, zz, cosmo)

        ll = ','.join([str(void_idx),str(ra), str(dec), str(z), str(this_reff)])

        gg.write(ll+'\n')


        void_idx += 1



mocks_file = '../mocks_all_info_w_smoothened.p'

mocks_data = pickle.load(open(mocks_file, 'rb'))

mocks_data = mocks_data['mocks_HR4_PH']

mocks_scale = 0.72
H0 = 72
Om0 = 0.26
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

for mock in tqdm(mocks_data):


    void_idx = 1

    gg.write(f'#{mock}\n')
    gg.write(f'#void_idx,ra,dec,z,Reff\n')

    for this_boundary in mocks_data[mock]['smoothened_minimal']:

        this_reff = hf.get_reff_convexhull(this_boundary)*mocks_scale

        centroid = np.mean(this_boundary, axis=0)
        xx = centroid[0]
        yy = centroid[1]
        zz = centroid[2]
        ra, dec, z = hf.invert(xx, yy, zz, cosmo)

        ll = ','.join([str(void_idx),str(ra), str(dec), str(z), str(this_reff)])

        gg.write(ll+'\n')

        void_idx += 1


gg.close()


