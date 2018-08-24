from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from explorator.imports.maps.map_standard4d_v3 import Mapping

def main():

    mapping = Mapping(k1 = -0.5,k2 = -0.5, k =0.01)
    num_iter = 10**5
    point = np.random.rand(4)
    #point = np.array([ 0.99212823,  .99212823,  0.81481635,  0.51498974])
    orbit = mapping.mapN(point, num_iter)
    x1 = orbit.points[:,0]
    x2 = orbit.points[:,1]
    x3 = orbit.points[:,2]
    x4 = orbit.points[:,3]

    plt.plot(x1-x2)
    plt.show()


###############################################################################
# Call `main()` if run from the command line
###############################################################################

if __name__ == "__main__":
    main()