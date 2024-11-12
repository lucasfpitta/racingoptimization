import numpy as np
import matplotlib.pyplot as plt


#Force matrix R translated to the path - R_t. 
#Input angle at midpoints, wheels angles,vehicle info
def force_tilde(angles,theta_r,theta_f0,theta_f1,n_wheels,width,L,Wf,h):
    if n_wheels!=3:
        print("incompatible number of wheels, check force_tilde at model4")
        SystemExit
        
    
    R_t = np.zeros((len(angles),6,9))
    
    #iteration over each section
    for i in range(len(angles)):   
        R_t[i]= np.array([[\
            #first coordinate
            np.cos(angles[i]+theta_f0[i]), -np.sin(angles[i]+theta_f0[i]),0,\
            np.cos(angles[i]+theta_f1[i]),-np.sin(angles[i]+theta_f1[i]),0,\
            np.cos(angles[i]-theta_r[i]),-np.sin(angles[i]-theta_r[i]),0], \
            #second coordinate        
            [np.sin(angles[i]+theta_f0[i]), np.cos(angles[i]+theta_f0[i]),0,\
            np.sin(angles[i]+theta_f1[i]),np.cos(angles[i]+theta_f1[i]),0,\
            np.sin(angles[i]-theta_r[i]),np.cos(angles[i]-theta_r[i]),0],\
             #third coordinate 
             [0,0,1,0,0,1,0,0,1],\
             #First angle
             [-width/2*np.cos(theta_f0[i]+theta_r[i])+Wf*L*np.sin(theta_f0[i]+theta_r[i]),\
              Wf*L*np.cos(theta_f0[i]+theta_r[i])+width/2*np.sin(theta_f0[i]+theta_r[i]),0,\
             width/2*np.cos(theta_f1[i]+theta_r[i])+Wf*L*np.sin(theta_f1[i]+theta_r[i]),
             Wf*L*np.cos(theta_f1[i]+theta_r[i])-width/2*np.sin(theta_f1[i]+theta_r[i]),0,\
                 0,-(1-Wf)*L,0],\
             #Second angle
             [-h*np.cos(theta_f0[i]+theta_r[i]),h*np.sin(theta_f0[i]+theta_r[i]),-Wf*L,\
              -h*np.cos(theta_f1[i]+theta_r[i]),h*np.sin(theta_f1[i]+theta_r[i]),-Wf*L,\
                  -h,0,(1-Wf)*L],
             #Third angle           
             [h*np.sin(theta_f0[i]+theta_r[i]),h*np.cos(theta_f0[i]+theta_r[i]),width/2,\
              h*np.sin(theta_f1[i]+theta_r[i]),h*np.cos(theta_f1[i]+theta_r[i]),-width/2,\
                  0,h,0]])
    return R_t












#Mass matrix M translated to the path - M_t. 
# Input is a scipy path derivative, angle derivative, a midpoint 
# discretization vector over [0,1], and vehicle info
def mass_tilde(derivative, angle_derivative, discretization,m,J):
    zero_lines = np.zeros(len(discretization))
    M_t = np.transpose(np.vstack((m*derivative(discretization),\
        zero_lines,J*angle_derivative,zero_lines,zero_lines)))
    return M_t








#Centrifugal matrix C translated to the path - C_t. 
# Input is a scipy path derivative and second derivative, angle second derivative,
#a midpoint discretization vector over [0,1], vehicle info
def centrifugal_tilde(derivative,secondderivative,angle_sec_derivative,\
    discretization,m,J,pho_air,A0,Cx):
    zero_lines = np.zeros(len(discretization))
    C_t = np.hstack((m*np.transpose(secondderivative(discretization))-pho_air*A0*Cx/2*\
        np.transpose(derivative(discretization))*(np.linalg.norm(np.transpose(\
            derivative(discretization)),axis=1)[:, np.newaxis]),
        zero_lines.reshape(-1, 1),J*angle_sec_derivative.reshape(-1, 1),\
            zero_lines.reshape(-1, 1),zero_lines.reshape(-1, 1)))
    return C_t








#Independet forces vector d 
# Input midpoint discretization vector over [0,1], vehicle info

def independent_tilde(discretization,m):
    d = np.array([0,0,m*9.81,0,0,0])
    return np.tile(d,(len(discretization),1))












#Power matrix A translated to the path - A_t. 
# Input is a scipy path derivative, a midpoint discretization 
# vector over [0,1], the number of wheels
def power_tilde(derivative,angles, discretization,theta_r,theta_f0,theta_f1):
    
    Rotation = np.zeros((len(angles),9,2))
    
    for i in range(len(angles)):
        Rotation[i]= np.array([[np.cos(angles[i]+theta_f0[i]), \
            np.sin(angles[i]+theta_f0[i])],\
            [-np.sin(angles[i]+theta_f0[i]),np.cos(angles[i]+theta_f0[i])],[0,0],\
            [np.cos(angles[i]+theta_f1[i]),np.sin(angles[i]+theta_f1[i])],\
            [-np.sin(angles[i]+theta_f1[i]),np.cos(angles[i]+theta_f1[i])],[0,0],\
            [np.cos(angles[i]-theta_r[i]),np.sin(angles[i]-theta_r[i])],\
            [-np.sin(angles[i]-theta_r[i]),np.cos(angles[i]-theta_r[i])],[0,0]])
    A_t = (Rotation@np.transpose(derivative(discretization))\
        [..., np.newaxis]).squeeze(-1)
    return A_t





#Defines the necessary vectors and matrix to model1
#Input scipy spline
#Output matrices
def model4(spline,angles,angle_derivative,angle_sec_derivative,\
    theta_r,theta_f0,theta_f1,\
    M,m,mu,pho_air,A0,Cx,J,width,L,Wf,h,n_wheels):

    
    #discretization lenght
    deltatheta = 1/(M-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = M-1)
    
    
    
    #Force matrix R translated to the path - R_t.
    R_t = force_tilde(angles,theta_r,theta_f0,theta_f1,n_wheels,width,L,Wf,h)
    
    #Other matrices
    M_t = mass_tilde(spline.derivative(),angle_derivative,discretization,m,J)
    C_t = centrifugal_tilde(spline.derivative(),spline.derivative().derivative(),\
        angle_sec_derivative,discretization, m,J,pho_air,A0,Cx)
    d_t = independent_tilde(discretization,m)
    A_t = power_tilde(spline.derivative(),angles,discretization,theta_r,theta_f0,theta_f1)
    return R_t, M_t, C_t, d_t, A_t