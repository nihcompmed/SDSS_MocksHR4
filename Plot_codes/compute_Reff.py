import numpy as np
from scipy.spatial import ConvexHull
import math


def get_reff(boundary):

    vol = ConvexHull(boundary).volume
    return eff_rad = (3*vol/(4*math.pi))**(1/3)





