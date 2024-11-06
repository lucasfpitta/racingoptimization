import numpy as np
import matplotlib.pyplot as plt


#Force matrix R translated to the path - R_t. 
#Input angle at midpoints
def force_tilde(angles):
    R_t= np.array([[[np.cos(theta), -np.sin(theta)], \
        [np.sin(theta), np.cos(theta)]] for theta in angles])
    return R_t












#Mass matrix M translated to the path - M_t. 
# Input is a scipy path derivative and a midpoint discretization 
# vector over [0,1]
def mass_tilde(derivative, discretization):
    M_t = np.transpose(derivative(discretization))
    return M_t








#Centrifugal matrix C translated to the path - C_t. 
# Input is a scipy path derivative and second derivative and 
#a midpoint discretization vector over [0,1]
def centrifugal_tilde(derivative,secondderivative, discretization,\
        pho_air,A0,Cx):
    C_t = np.transpose(secondderivative(discretization))-pho_air*A0*Cx/2*\
        np.transpose(derivative(discretization))**2
    return C_t







#Power matrix A translated to the path - A_t. 
# Input is a scipy path derivative and a midpoint discretization 
# vector over [0,1]
def power_tilde(derivative,angles, discretization):
    Rotation= np.array([[[np.cos(theta), -np.sin(theta)], \
        [np.sin(theta), np.cos(theta)]] for theta in angles])
    A_t = (Rotation@np.transpose(derivative(discretization))\
        [..., np.newaxis]).squeeze(-1)
    return A_t





#Defines the necessary vectors and matrix to model1
#Input scipy spline
#Output matrices
def model2(spline,angles,M):
    
    #discretization lenght
    deltatheta = 1/(M-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = M-1)
    
    
    m = 85 #vehicle mass
    mu = 1 #tyre friction coeficient 
    pho_air = 1.225 #air density
    A0 = 0.1 #frontal area of the car
    Cx = 0.1 #Drag coeficient
    
    
    #Force matrix R translated to the path - R_t.
    R_t = force_tilde(angles)
    
    #Other matrices
    M_t = m*mass_tilde(spline.derivative(), discretization)
    C_t = m*centrifugal_tilde(spline.derivative(),spline.derivative().derivative(),discretization,\
        pho_air,A0,Cx)
    A_t = power_tilde(spline.derivative(),angles,discretization)
    
    return R_t, M_t, C_t, A_t