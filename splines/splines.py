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
     
    #connects the last and first points    
    splpoints[0][-1]=splpoints[0][0]
    splpoints[1][-1]=splpoints[1][0]
    return splpoints








#constructs the scipy spline and scipy spline derivatives
#Input 2d control point vecor splpoints
#Output scipy spline interpol, and scipy spline derivative, diff

def splines_and_derivatives(splpoints):
    
    t = np.linspace(0,1,num = len(splpoints[0]))
    interpol = scp.interpolate.CubicSpline(t, (splpoints[0],splpoints[1]
                                            ),axis=1, bc_type='periodic')
    
    #spline derivative
    diff = interpol.derivative()

    return interpol, diff








#define angle assessments over the trajectory. 
#Input assessments of the derivative of the spline in respect to theta 
# over the trajecory 
#Output angle assessment vector


def find_angle(spline,N_angle):
    
    #defines spline info
    delta  = 1/(N_angle-1)
    midpoints = np.linspace(delta/2,1-delta/2,num = \
                    (N_angle-1))
    
    #spline derivatives
    diff = spline.derivative()
    second_diff = spline.derivative().derivative()
    third_diff = spline.derivative().derivative().derivative()
    
    #splines at the midpoints
    spline_points = spline(midpoints)
    derivative = diff(midpoints)
    sec_derivative=second_diff(midpoints)
    third_derivative=third_diff(midpoints)
    
    
    #Builds the vectors
    angle = np.zeros(len(midpoints))
    angle_derivative = np.zeros(len(midpoints))
    angle_sec_derivative = np.zeros(len(midpoints))
    
    
    #calculates the angle with the spline derivative
    for i in range(len(midpoints)):
        if np.abs(derivative[0][i]) <= 1e-18 and np.abs(derivative[1][i]) <= 1e-18:
            raise ValueError("Valores de x` e y` estão muito próximos de zero!")
        elif np.abs(derivative[0][i]) <= 1e-18:
            angle[i] = np.pi/2
        elif derivative[0][i] >= 0:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))
        else:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))+np.pi
            
            
    #calculates the  derivative
    angle_derivative = (derivative[0]*sec_derivative[1]-derivative[1]*sec_derivative[0])\
        /(derivative[0]**2+derivative[1]**2)
    
    
    #calculates second derivative (complicated expression)
    angle_sec_derivative = 1/(derivative[0]**2+derivative[1]**2)**2*\
        ((derivative[0]**2+derivative[1]**2)*(derivative[0]*\
            third_derivative[1]-derivative[1]*third_derivative[0])-\
            (derivative[0]*sec_derivative[1]-derivative[1]*sec_derivative[0])\
        *(2*derivative[0]*sec_derivative[0]+2*derivative[1]*sec_derivative[1]))   
    
    return angle, angle_derivative, angle_sec_derivative



#defines the front and rear wheels angles in the trajectory
def model4_extra_angles(spline_derivative,spline_sec_derivative,\
    n_discretization,Wf,L,w):
    
    #midpoints
    delta = 1/(n_discretization-1)
    discretization = np.linspace(delta/2,1-delta/2,num=n_discretization-1)
    
    #Cg instantaneous radius
    Rc = (spline_derivative(discretization)[0]**2+spline_derivative(\
        discretization)[1]**2)**(3/2)/(spline_derivative(discretization)[0]*\
        spline_sec_derivative(discretization)[1]-spline_sec_derivative(\
            discretization)[0]*spline_derivative(discretization)[1])
        
    #Rear axel instantaneus radious    
    Rr = Rc/np.abs(Rc)*np.sqrt(Rc**2-((1-Wf)*L*np.ones(len(discretization))\
        )**2)
    
    #rear wheel angle
    theta_r = np.arcsin((1-Wf)*L/Rc)
    
    #front wheel angles
    theta_f0 = np.arctan(L/(Rr-w/2))-theta_r
    theta_f1 = np.arctan(L/(Rr+w/2))-theta_r
    
    return theta_r,theta_f0,theta_f1
    






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
    
    angle, angle_derivative, angle_sec_derivative  = \
        find_angle(spline,N_angle)
    
    return spline, derivative, angle, angle_derivative, angle_sec_derivative