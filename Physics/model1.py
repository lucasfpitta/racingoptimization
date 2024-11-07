import numpy as np



#Mass matrix M translated to the path - M_t. 
# Input is a scipy path derivative and a midpoint discretization 
# vector over [0,1]
def mass_tilde(derivative, discretization):
    M_t = np.transpose(derivative(discretization))
    return M_t








#Centrifugal matrix C translated to the path - C_t. 
# Input is a scipy path second derivative and a midpoint discretization
# vector over [0,1]
def centrifugal_tilde(secondderivative, discretization):
    C_t = np.transpose(secondderivative(discretization))
    return C_t







#Power matrix A translated to the path - A_t. 
# Input is a scipy path derivative and a midpoint discretization 
# vector over [0,1]
def power_tilde(derivative, discretization):
    A_t = np.transpose(derivative(discretization))
    return A_t





#Defines the necessary vectors and matrix to model1
#Input scipy spline
#Output matrices
def model1(spline,M,m,mu):
    
    #discretization lenght
    deltatheta = 1/(M-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = M-1)
    
    
    #Force matrix R translated to the path - R_t.
    R_t = np.tile(np.identity(2),(M-1,1,1))
    
    
    #Other matrices
    M_t = m*mass_tilde(spline.derivative(), discretization)
    C_t = m*centrifugal_tilde(spline.derivative().derivative(),discretization)
    A_t = power_tilde(spline.derivative(), discretization)
    
    return R_t, M_t, C_t, A_t