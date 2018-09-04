"""The Smaller Alignment Index (SALI) computation.

References:
Sko2001
SkoMan2016
"""

from __future__ import division, print_function

import numpy as np


##############################################################################
# SALI Computation
##############################################################################

def compute_sali(initial, renorm_time, max_time):
    """SALI computation.

    Parameters
    ----------
    initial - initial point
    mapping - the given map
    renorm_time - the renormalization time
    max_time - the maxumum number of iterations

    Returns
    -------
    sali - value of computed SALIs
    n_time - iterations where SALIs are computed

    """
    dim = 2
    # Threshold value for stopping the computation.
    # If the SALI goes below this, the trajectory is considered
    # to be chaotic.
    # Note: there are examples in the literatur FIXME Ref,
    # showing that after some decrease the SALI could become larger again.
    # FIXME: this should be a parameter of the function.
    sali_thrsld = 10**-16
    print("Defining the threshold to be ", sali_thrsld)
    grm = np.random.rand(dim, 2)  # random deviation vector

    # Orthonormalization using Gram-Schmidt procedure
    w1_initial = grm[:, 0]

    w2_initial = grm[:, 1] - (
        np.dot(grm[:, 1], w1_initial)*w1_initial /
        np.linalg.norm(w1_initial)**2)

    w1_initial = w1_initial / np.linalg.norm(w1_initial)
    w2_initial = w2_initial / np.linalg.norm(w2_initial)

    # Initial deviation matrix.
    # w1 and w1 are initial deviation vectors and are orthogonal
    # FIXME AB: should this be ``w1 and w2 ...''?
    w1_n = w1_initial
    w2_n = w2_initial
    sali1 = np.linalg.norm(w1_n + w2_n)
    sali2 = np.linalg.norm(w1_n - w2_n)
    sali_diff = min(sali1, sali2)
    sali = []      # computed SALI is stored here
    sali.append(sali_diff)
    n_time = []
    n_time.append(1)

    for i in range(max_time):
        jacobian = Jacob(initial)
        w1_n = np.dot(jacobian, w1_n)
        w2_n = np.dot(jacobian, w2_n)

        if i % renorm_time == 0:
            w1_n = w1_n / np.linalg.norm(w1_n)
            w2_n = w2_n / np.linalg.norm(w2_n)
            sali1 = np.linalg.norm(w1_n + w2_n)
            sali2 = np.linalg.norm(w1_n - w2_n)
            sali_diff = min(sali1, sali2)
            sali.append(sali_diff)
            n_time.append(i+1)  # i starts from zero

            if sali_diff < sali_thrsld:
                break
        initial = sna_map(initial)



    return sali, n_time

##############################################################################

def sna_map(initial):


    sigma = 1.2
    omega = 0.618
    x = initial[0]
    theta = initial[1]

    x1 = 2.0*sigma*np.tanh(x)*np.cos(2.0*np.pi*theta)
    theta1 = np.mod(theta + omega,1)

    initial1 = [x1,theta1]

    return initial1

def Jacob(initial, sigma = 1.2):

    x = initial[0]
    theta = initial[1]

    J = np.zeros([2,2])
    # from IPython import embed
    # embed()

    J[0,0] = 2.0*sigma*(1.0 - np.tanh(x)**2)*np.cos(2.0*np.pi*theta)
    J[0,1] = 2.0*sigma*np.tanh(x)*(-2.0*np.pi)*np.sin(2.0*np.pi*theta)

    J[1,0] = 0.0
    J[1,1] = 1.0

    return (J)

def main():

    initial = np.random.rand(2)

    sali, ntime = compute_sali(initial=initial, renorm_time= 10,
                               max_time=10**5)
    print(sali)
    print(ntime)

# Call `main()` if run from the command line
###############################################################################
if __name__ == "__main__":
    main()
###############################################################################