from __future__ import division, print_function
import numpy as np
import matplotlib.animation as animation
from scipy.integrate import odeint
from numpy import arange
from scipy.linalg import block_diag
from CPG.hdf import HDFClassIO
import CPG
from pylab import *
from sympy.physics.quantum.tensorproduct import TensorProduct
def coherent_std_map(eps = 0.0, delta = 0.1, K = 6.0):

    N = 10**7
    n_incod = 50
    del_Q = np.zeros(N)
    P1_incod = np.random.uniform(-0.5,0.5,n_incod)
    Q1_incod = np.random.uniform(0,1, n_incod)
    Q2_incod = np.random.uniform(0,1,n_incod)

    P_d = []
    Q_d = []
    Q_r = []
    sync_time_control = []
    for j in range(n_incod):
        P1 = np.zeros(N)
        P2 = np.zeros(N)
        Q1 = np.zeros(N)
        Q2 = np.zeros(N)
        P1[0] = P1_incod[j]
        Q1[0] = Q1_incod[j]
        Q2[0] = Q2_incod[j]
        for i in range(N-1):

            if abs(P1[i]) < delta and abs(Q1[i]-0.5) < delta:
                K1 = K - eps
            else:
                K1 = K
            #print(K1)
            P1[i+1] = P1[i] + K1/(2.0*np.pi)*np.sin(2.0*np.pi*(Q1[i]))
            # P1[i+1] = (P1[i+1] + 0.5) % 1 - 0.5
            Q1[i+1] = P1[i+1] + Q1[i]

            if abs(P1[i]) < delta and abs(Q2[i]-0.5) < delta:
                K2 = K - eps
            else:
                K2 = K

            P2[i+1] = P1[i] + K2/(2.0*np.pi)*np.sin(2.0*np.pi*(Q2[i]))
            Q2[i+1] = P2[i+1] + Q2[i]

            P1[i+1] = (P1[i+1] + 0.5) % 1 - 0.5
            Q1[i+1] = Q1[i+1] % 1
            Q2[i+1] = Q2[i+1] % 1

        diff = np.abs(Q1 - Q2)
        # from IPython import embed
        # embed()

        for i in range(N-10):
            if diff[i] < 10**-5:
                # from IPython import embed
                # embed()
                test = diff[i+1:N]
                if test.max() < 10**-5:
                    sync_time_control.append(i)
                    #print(i)
                    break

        P_d.append(P1)
        Q_d.append(Q1)
        Q_r.append(Q2)

    print(sync_time_control)

    # without control

    P_d = []
    Q_d = []
    Q_r = []
    sync_time = []
    for j in range(n_incod):
        P1 = np.zeros(N)
        P2 = np.zeros(N)
        Q1 = np.zeros(N)
        Q2 = np.zeros(N)
        P1[0] = P1_incod[j]
        Q1[0] = Q1_incod[j]
        Q2[0] = Q2_incod[j]
        for i in range(N - 1):
            P1[i + 1] = P1[i] + K1 / (2.0 * np.pi) * np.sin(
                2.0 * np.pi * (Q1[i]))
            # P1[i+1] = (P1[i+1] + 0.5) % 1 - 0.5
            Q1[i + 1] = P1[i + 1] + Q1[i]
            P2[i + 1] = P1[i] + K1 / (2.0 * np.pi) * np.sin(
                2.0 * np.pi * (Q2[i]))
            Q2[i + 1] = P2[i + 1] + Q2[i]

            P1[i + 1] = (P1[i + 1] + 0.5) % 1 - 0.5
            Q1[i + 1] = Q1[i + 1] % 1
            Q2[i + 1] = Q2[i + 1] % 1

        diff = np.abs(Q1 - Q2)
        # from IPython import embed
        # embed()

        for i in range(N - 10):
            if diff[i] < 10 ** -5:
                # from IPython import embed
                # embed()
                test = diff[i + 1:N]
                if test.max() < 10 ** -5:
                    sync_time.append(i)
                    #print(i)
                    break

    #diff_Q = abs(np.array(Q_d) - np.array(Q_r))
    # from IPython import embed
    # embed()
    # plt.plot(Q_d,P_d,'.b')
    # #plt.plot(diff_Q[0],'b')
    # plt.show()

    print(sync_time)

    data = dict(without_control=np.array(sync_time), with_control=np.array(sync_time_control))
    fname = "sync_times.hf"
    CPG.hdf.save_h5_dict(fname, data)

if __name__ == "__main__":
    coherent_std_map(eps = 2, delta = 0.2, K = 6.0)
    #test()


