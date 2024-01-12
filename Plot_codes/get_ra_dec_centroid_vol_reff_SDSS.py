from astropy.cosmology import FlatLambdaCDM
import astropy.cosmology.units as cu
import astropy.units as u
import astropy
import math
import numpy as np

#print(astropy.__version__)
cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
#print(cosmo)

def polar_to_cart(a, d, z):

    # Right ascension alpha, a
    # Declination, d
    # distance z

    # Convert a and d from degrees to radians
    a = a * math.pi/180
    d = d * math.pi/180

    # Cartesian (xx, yy, zz)

    xx = z * np.cos(a) * np.cos(d)
    yy = z * np.sin(a) * np.cos(d)
    zz = z * np.sin(d)

    return xx, yy, zz

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






    ##print(math.atan(yy/xx))
    #print(math.asin(yy/math.sqrt(val)))
    #print(math.acos(xx/math.sqrt(val)))


    ## Right ascension alpha, a
    ## Declination, d
    ## distance z
    #
    ## Convert a and d from degrees to radians
    #a = a * math.pi/180
    #d = d * math.pi/180
    #
    ## Cartesian (xx, yy, zz)
    #
    #xx = z * np.cos(a) * np.cos(d)
    #yy = z * np.sin(a) * np.cos(d)
    #zz = z * np.sin(d)

    return a, d, z


#exit()
# File 1: Raw original data
fname = '../V2_analysis/input_galaxy_Mr20.19_z0.02-0.116_redshift'
ff = open(fname+'.cat', 'r')

max_diff_ra = 0
max_diff_dec = 0
max_diff_z = 0

count = 0

for line in ff:


    line = line.strip('\n')
    line = line.split(' ')
    ra = float(line[0])
    dec = float(line[1])
    z = float(line[2])

    weight = float(line[3])

    if weight == 0:
        continue


    z = z * cu.redshift

    z_comoving = cosmo.comoving_distance(z)
    #print(z_comoving)
    #exit()

    #d = z.to(u.Mpc, cu.redshift_distance("FlatLambdaCDM"\
    #                                    , kind="comoving"\
    #                                    , H0=100\
    #                                    , Om0=0.3\
    #                                    ))

    #print(z)
    #print(d)

   #print(ra, dec, z)

    if z_comoving.value < 0:
        print(z_comoving)
        exit()

    xx, yy, zz = polar_to_cart(ra, dec, z_comoving.value)

    ra_hat, dec_hat, z_hat = invert(xx, yy, zz)

    max_diff_ra = max(max_diff_ra, abs(ra-ra_hat))
    max_diff_dec = max(max_diff_dec, abs(dec-dec_hat))
    max_diff_z = max(max_diff_z, abs(z-z_hat))

    #max_diff = max(max_diff, abs(ra-ra_hat), abs(z-z_hat), abs(dec-dec_hat))

    print(count, max_diff_ra, max_diff_dec, max_diff_z, end='\r')
    count += 1


    #print(z)
exit()




ff.close()

                                   
