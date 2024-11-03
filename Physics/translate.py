import numpy as np




#translate generalized velocity b into cartesian velocity v
#Input derivative of the spline, b, and the number of path sections
def translate_velocity(spline_diff, b,N_discretization):
    
    #discretization lenght
    deltatheta = 1/(N_discretization-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,
                               num = N_discretization-1)
    spline_derivative = spline_diff(discretization)
    
    #create velocity v
    v = np.zeros(N_discretization-1)
    
    for i in range(N_discretization-1):
        v[i] = np.sqrt((spline_derivative[0][i]*np.sqrt((b[i]+b[i+1])/2)
                )**2+(spline_derivative[1][i]*np.sqrt((b[i]+b[i+1])/2))**2)
        
    return  v








#translate generalized acceleration into cartesian acceleration a
#Input derivative of the spline, second derivative, b, a, and the 
# number of path sections

def translate_acceleration(spline_diff, spline_second_diff, b, a,N_discretization):
    
    #discretization lenght
    deltatheta = 1/(N_discretization-1)
    
    #midpoints discretization
    discretization=np.linspace(deltatheta/2,1-deltatheta/2,num = N_discretization-1)
    
    #calculates derivatives at the points
    spline_derivative = spline_diff(discretization)
    spline_second_derivative = spline_second_diff(discretization)
    
    
    #create acceleration a
    a = np.zeros(N_discretization-1)
    for i in range(N_discretization-1):
        a[i] = np.sqrt((spline_derivative[0][i]*a[i])**2+(spline_derivative[1][i]*a[i]
                )**2)+ np.sqrt((spline_second_derivative[0][i]*(b[i]+b[i+1])/2)**2+(
                    spline_second_derivative[1][i]*(b[i]+b[i+1])/2)**2)
            
    return  a
