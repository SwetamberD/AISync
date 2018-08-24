from __future__ import division, print_function
import numpy as np
from numpy import arange
from CPG.hdf import HDFClassIO
import CPG
import math

def coherent_std_map(P1_incod, Q1_incod, Q2_incod, eps = 0.0, delta = 0.1, K = 6.0, N = 10, n_incod = 10):
    sync_time = []
    for j in range(n_incod):
        P1 = np.zeros(N)
        P2 = np.zeros(N)
        Q1 = np.zeros(N)
        Q2 = np.zeros(N)
        P1[0] = P1_incod[j]
        Q1[0] = Q1_incod[j]
        Q2[0] = Q2_incod[j]
        for i in range(N-1):
            if abs(P1[i]) < delta and abs(Q1[i]) < delta:
                K1 = K - eps
            else:
                K1 = K
            P1[i+1] = P1[i] + K1/(2.0*math.pi)*math.sin(2.0*math.pi*Q1[i])
            #P1[i+1] = (P1[i+1] + 0.5) % 1 - 0.5
            Q1[i+1] = P1[i+1] + Q1[i]

            if abs(P1[i]) < delta and abs(Q2[i]) < delta:
                K2 = K - eps
            else:
                K2 = K
            P2[i+1] = P1[i] + K2/(2.0*math.pi)*math.sin(2.0*math.pi*Q2[i])
            Q2[i+1] = P2[i+1] + Q2[i]

            if P1[i+1] > 0.5:
                P1[i+1] = P1[i+1] - 1.0
            elif P1[i+1] < -0.5:
                P1[i+1] = P1[i+1] + 1.0
            Q1[i+1] = Q1[i+1] % 1.0
            Q2[i+1] = Q2[i+1] % 1.0
            if i % 1000 == 0 and i != 0:
               # print("i =",i)
                # from IPython import embed
                # embed()
                diff = np.abs(Q1[i-1000:i] - Q2[i-1000:i])
                if diff.max() < 10 ** -5:
                    break
        # P_d.append(P1)
        # Q_d.append(Q1)
        # Q_r.append(Q2)
        # P_d.append(P1[N - 10000:N])
        # Q_d.append(Q1[N - 10000:N])
        # Q_r.append(Q2[N - 10000:N])
        diff = np.abs(Q1 - Q2)
        # from IPython import embed
        # embed()
        t = 0
        len_diff = len(diff)
        for i in range(len_diff):
            if diff[i] < 10**-5:
                # from IPython import embed
                # embed()
                test = diff[i+1:-1]
                if test.max() < 10**-5:
                    t = i
                    sync_time.append(t)
                    print("sync time = ", i)

                    break
        if t == 0:
            sync_time.append(t)
            print("sync time = ", t)

    # plt.plot(Q_d,P_d,'.b')
    # plt.xlim(0,1)
    # plt.ylim(-0.5,0.5)
    # plt.show()
    #print(sync_time_control)
    return(sync_time)


    # data = dict(without_control=np.array(sync_time), with_control=np.array(sync_time_control))
    # fname = "sync_times_in_cod_50000.hf"
    # CPG.hdf.save_h5_dict(fname, data)

def sync_code(save_data = False):
    N = 10 ** 7
    n_incod = 10
    del_Q = np.zeros(N)
    P1_incod = np.random.uniform(-0.4, 0.4, n_incod)
    Q1_incod = np.random.uniform(0.3, 0.7, n_incod)
    Q2_incod = np.random.uniform(0.3, 0.7, n_incod)

    sync_time_without_control = coherent_std_map(N = N, n_incod = n_incod, eps = 0.0, delta = 0.0, K = 6, P1_incod = P1_incod, Q1_incod = Q1_incod, Q2_incod = Q2_incod)
    sync_time_with_control = coherent_std_map(N = N, n_incod = n_incod, eps = 1.5, delta = 0.4, K = 6, P1_incod = P1_incod, Q1_incod = Q1_incod, Q2_incod = Q2_incod)

    print(sync_time_without_control)
    print(sync_time_with_control)
    if save_data:
        data = dict(without_control=np.array(sync_time_without_control), with_control=np.array(sync_time_with_control))
        fname = "sync_times_in_cod_50000.hf"
        CPG.hdf.save_h5_dict(fname, data)

def do_plots():
    fname = "sync_times.hf"
    data = CPG.hdf.load_h5_dict(fname)

    time_c = data["with_control"]
    time_w = data["without_control"]

    from IPython import embed
    embed()



if __name__ == "__main__":
    #coherent_std_map(eps = 1.5, delta = 0.4, K = 6)
    sync_code()
    #test()
    # do_plots()


