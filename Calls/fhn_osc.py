from __future__ import division, print_function
import numpy as np
import matplotlib.animation as animation
from scipy.integrate import odeint
from numpy import arange
from scipy.linalg import block_diag

from pylab import *
from sympy.physics.quantum.tensorproduct import TensorProduct
def fitzhugh_nagumo(state,t, L = -0.2, eps = 0.1, D = 0.1):
    u,w = state
    print(L)
    du = -w -u*(u-1)*(u-L)
    dw = eps*(u-D*w)
    return [du, dw]


def rossler(state,t, a = 0.2, b = 0.2, c = 9.0 ):
    x,y,z = state

    x_dot = -y - z
    y_dot = x + a*y
    z_dot = b + (x-c)*z

    return [x_dot, y_dot, z_dot]
def main1():
    t = arange(0, 60, 0.01)
    init_state = [.1, .1]
    L = -0.5
    eps= 0.1
    D = 0.1
    state = odeint(fitzhugh_nagumo, init_state, t, args=(L,eps,D))
    fig = figure()
    xlabel('u')
    ylabel('w')
    plot(state[:, 0], state[:, 1], 'b-', alpha=0.8)
    plt.show()

def test_finest_simul_block_diag():
    # Given matrices
    A  = np.array([[2.0,1.0, 0.0,0.0],
                 [1.0,2.0,0.0,0.0],
                 [0.0,0.0,1.0,2.0],
                 [0.0,0.0,2.0,1.0]])

    B = np.array([[0.0,0.0, 1.0,0.0],
                 [0.0,0.0,0.0,1.0],
                 [1.0,0.0,0.0,0.0],
                 [0.0,1.0,0.0,0.0]])

    C = np.array([[0.0,0.0, 0.0,1.0],
                 [0.0,0.0,1.0,0.0],
                 [0.0,1.0,0.0,0.0],
                 [1.0,0.0,0.0,0.0]])

    # identity matrix
    id = np.eye(4)

    # constructing n^2 by n^2 matrices for A, B, and C using id
    T1 = TensorProduct(id, A) - TensorProduct(A.transpose(), id)
    T2 = TensorProduct(id, B) - TensorProduct(B.transpose(), id)
    T3 = TensorProduct(id, C) - TensorProduct(C.transpose(), id)

    # constructing the S matrix (see the algorithm)
    des_sum = np.matmul(T1.transpose(),T1) + np.matmul(T1, T1.transpose()) \
                + np.matmul(T2.transpose(),T2) + np.matmul(T2, T2.transpose()) \
                + np.matmul(T3.transpose(),T3) + np.matmul(T3, T3.transpose())

    # P = np.array([[-0.5, -0.5, -0.433, -0.559],
    #        [-0.5, -0.5, 0.433, 0.559],
    #        [-0.5, 0.5, -0.559, 0.433],
    #        [-0.5, 0.5, 0.559, -0.433]])

    eigenValues, eigenVectors = linalg.eig(des_sum)


    # sorting the eigenvalues and eigenvectors to find smallest eigs
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]

    # separating eig vectors corresponding to smallest eig values
    vec_A = eigenVectors[:,13]
    vec_B = eigenVectors[:,14]
    vec_C = eigenVectors[:,15]

    # constructing random coefficients for the matrix u
    coeffs_1= np.random.rand(3)
    coeffs_norm = coeffs_1/sum(coeffs_1)
    coeffs = sqrt(coeffs_norm)

    vec_u = coeffs[0]*vec_A + coeffs[1]*vec_B+ coeffs[2]*vec_C

    mat_X = np.zeros([4,4])
    mat_X[:,0] = vec_u[0:4]
    mat_X[:,1] = vec_u[4:8]
    mat_X[:,2] = vec_u[8:12]
    mat_X[:,3] = vec_u[12:16]

    # Hermitial matrix form mat X
    Herm_X = 0.5*(mat_X+mat_X.transpose())
    eigvals, eigvecs = np.linalg.eigh(Herm_X)

    # Block Digonalizing A, B, and C

    Diag_A = np.matmul(eigvecs.transpose(),np.matmul(A,eigvecs))
    Diag_B = np.matmul(eigvecs.transpose(),np.matmul(B,eigvecs))
    Diag_C = np.matmul(eigvecs.transpose(),np.matmul(C,eigvecs))

    # Print the results
    print("Matrix X:")
    print(Herm_X)
    print("\n Unitary matrix P (to be used for block diagonalizing A,B, and C):")
    print(np.matrix.round(eigvecs))
    print("\n Block diagonal A:")
    print(np.matrix.round(Diag_A,4))
    print("\n Block diagonal B:")
    print(np.matrix.round(Diag_B,4))
    print("\n Block diagonal C:")
    print(np.matrix.round(Diag_C,4))
    from IPython import embed
    embed()


if __name__ == "__main__":
    test_finest_simul_block_diag()
    #test()


