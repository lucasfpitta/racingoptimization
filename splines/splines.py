import numpy as np
import scipy as scp






#Define the control points over the trajectory. 
#Input 2d outlines left and right, 1d weight vector alfas
#Output 2d control points vector splpoints

def knot_points(left, right, alfas):
    
    splpoints = np.zeros((2,len(right[0])))
    
    for i in range(len(right[0])):
        splpoints[0][i] = right[0][i]+alfas[i]*(left[0][i]-right[0][i])
        splpoints[1][i] = right[1][i]+alfas[i]*(left[1][i]-right[1][i])
    return splpoints








#constructs the scipy spline and scipy spline derivatives
#Input 2d control point vecor splpoints
#Output scipy spline interpol, and scipy spline derivative, diff

def splines_and_derivatives(splpoints):
    
    t = np.linspace(0,1,num = len(splpoints[0]))
    
    interpol = scp.interpolate.CubicSpline(t, (splpoints[0],splpoints[1]
                                            ),axis=1, bc_type='natural')
    
    #spline derivative
    diff = interpol.derivative()

    return interpol, diff








#define angle assessments over the trajectory. 
#Input assessments of the derivative of the spline in respect to theta 
# over the trajecory 
#Output angle assessment vector


def find_angle(derivative):
    
    angle = np.zeros(len(derivative[0]))
    
    #calculates the angle with the spline derivative
    for i in range(len(derivative[0])):
        if derivative[0][i] >= 0:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))
        else:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))+np.pi
            
     #last angle is the same as the first       
    angle[0]=angle[-1]
    
    return angle








#Defines path data
#Inputs, outline 2d vectors left and right, 1d alfas weight vector, number of 
# angle assessments N_angle
#Output scipy spline and derivative, angle assessments

def path_info(left, right, alfas,N_angle):
    
    #calculate controlpoints
    splpoints = knot_points(left, right, alfas)
    
    #fit a spline on the control points and calculate derivative
    spline, derivative = splines_and_derivatives(splpoints)
    
    #calculates spline angles at midpoints
    delta  = 1/(N_angle-1)
    angle = find_angle(derivative(np.linspace(delta/2,1-delta/2,num = \
                    (N_angle-1))))
    
    return spline, derivative, angle