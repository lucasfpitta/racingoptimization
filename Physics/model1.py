import numpy as np


def mass_tilde(derivative, discretization):
    M_t = np.transpose(derivative(discretization))
    return M_t
def centrifugal_tilde(secondderivative, discretization):
    C_t = np.transpose(secondderivative(discretization))
    return C_t

def power_tilde(derivative, discretization):
    A_t = np.transpose(derivative(discretization))
    return A_t

def model1(spline,M):
    deltatheta = 1/(M-1)
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = M-1)
    m = 30
    mu = 1
    R_t = np.tile(np.identity(2),(M-1,1,1))
    M_t = m*mass_tilde(spline.derivative(), discretization)
    C_t = m*centrifugal_tilde(spline.derivative().derivative(),discretization)
    A_t = power_tilde(spline.derivative(), discretization)
    return R_t, M_t, C_t, A_t