import numpy as np
import matplotlib.pyplot as plt


#Force matrix R translated to the path - R_t. 
#Input angle at midpoints
def force_tilde(angles,n_wheels,width,L,Wf):
    if n_wheels!=4:
        print("incompatible number of wheels, check force_tilde at model3")
        SystemExit
        
    #torque line of the matrix    
    torque = np.array([-width/2, Wf*L, width/2, Wf*L,\
        -width/2,-(1-Wf)*L, width/2, -(1-Wf)*L])
    
    R_t= np.array([[[np.cos(theta), -np.sin(theta)]*n_wheels, \
        [np.sin(theta), np.cos(theta)]*n_wheels,torque] for theta in angles])
    return R_t












#Mass matrix M translated to the path - M_t. 
# Input is a scipy path derivative, angle derivative, a midpoint 
# discretization vector over [0,1], and vehicle info
def mass_tilde(derivative, angle_derivative, discretization,m,J):
    M_t = np.transpose(np.vstack((m*derivative(discretization),\
        J*angle_derivative)))
    return M_t








#Centrifugal matrix C translated to the path - C_t. 
# Input is a scipy path derivative and second derivative, angle second derivative,
#a midpoint discretization vector over [0,1], vehicle info
def centrifugal_tilde(derivative,secondderivative,angle_sec_derivative,\
    discretization,m,J,pho_air,A0,Cx):
    C_t = np.hstack((m*np.transpose(secondderivative(discretization))+pho_air*A0*Cx/2*\
        np.transpose(derivative(discretization))*(np.linalg.norm(np.transpose(\
            derivative(discretization)),axis=1)[:, np.newaxis]),\
                J*angle_sec_derivative.reshape(-1, 1)))
    return C_t







#Power matrix A translated to the path - A_t. 
# Input is a scipy path derivative, a midpoint discretization 
# vector over [0,1], the number of wheels
def power_tilde(derivative,angles, discretization,n_wheels):
    Rotation= np.array([np.vstack([np.array([[np.cos(theta), np.sin(theta)], \
        [-np.sin(theta), np.cos(theta)]])]*n_wheels) for theta in angles])
    A_t = (Rotation@np.transpose(derivative(discretization))\
        [..., np.newaxis]).squeeze(-1)
    return A_t





#Defines the necessary vectors and matrix to model1
#Input scipy spline
#Output matrices
def model3(spline,angles,angle_derivative,angle_sec_derivative,\
    M,m,mu,pho_air,A0,Cx,J,width,L,Wf,n_wheels):

    
    #discretization lenght
    deltatheta = 1/(M-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = M-1)
    
    
    
    #Force matrix R translated to the path - R_t.
    R_t = force_tilde(angles,n_wheels,width,L,Wf)
    
    #Other matrices
    M_t = mass_tilde(spline.derivative(),angle_derivative,discretization,m,J)
    C_t = centrifugal_tilde(spline.derivative(),spline.derivative().derivative(),\
        angle_sec_derivative,discretization, m,J,pho_air,A0,Cx)
    A_t = power_tilde(spline.derivative(),angles,discretization,n_wheels)

    return R_t, M_t, C_t, A_t