from Map_processing.mapreader import get_outline
from Map_processing.point_reduction_and_smooting import reduce_points, smooth_moving_average
from splines.splines import path_info
import numpy as np









#Path imported by google eath. external = external outline google_earth 
# .kml points, 
#Inputs internal = internal outline google_earth .kml points, N_angle = 
#numer of angle assesments over the trajectory 
#Output scipy spline and derivative, outlines 2d vectors left and right, 
# 1d vector angle assessments angle.


def google_earth_path(external,internal,N_angle):
    
    alfa = 0.01 #smoothing parameter for elevation
    n_neighborhood = 16 #smoothing parameter for elevation
    max_slope = 0.06 #smoothing parameter for elevation
    
    
    
    #read the map and set the metric coordinates
    right, left,  = get_outline(external,internal)
    
    
    #Uncomment if inner and outer lines have different number of points
    #coords= reduce_points(right,left) 
    
    #smooth array 2 and 5 (z coordinate)
    # coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],
    #                     int(len(coords[0])/n_neighborhood),alfa,max_slope)
    
    
    
    #create a random path on the track
    alfas = np.random.random_sample(len(right[0]))
    alfas[-1]=alfas[0]
    
    
    #calculate the spline (scipy spline), its derivative (scipy spline), 
    # angle (array over N_angle points)
    #left and right are the outlines and alfas is the random path
    spline, derivative, angle = path_info(left, right, alfas,N_angle)
    return spline, derivative, angle, right, left














#defines a circular path with the vehicle passing on the middle of the track
# N_angle = numer of angle assesments over the trajectory 
#Output scipy spline and derivative, outlines 2d vectors left and right, 
# 1d vector angle assessments angle.

def circle(N_angle):
    
    Radius=100
    Delta_radius=10
    
    #define track middle
    alfas= alfas = np.ones(30)*0.5
    theta = np.linspace(0, 2 * np.pi, len(alfas))
    
    
    #outer line
    right = [(Radius+Delta_radius)*np.cos(theta),(Radius+Delta_radius)
             *np.sin(theta)]
    
    #inner line
    left = [(Radius-Delta_radius)*np.cos(theta),(Radius-Delta_radius)
            *np.sin(theta)]
    
    #calculate the spline (scipy spline), its derivative (scipy spline), 
    # angle (array over N_angle points)
    #left and right are the outlines and alfas is the random path
    spline, derivative, angle = path_info(left, right, alfas,N_angle)
    
    
    return spline, derivative, angle, right, left












#defines a semi_circular path with the vehicle passing on the middle of the track
# N_angle = numer of angle assesments over the trajectory 
#Output scipy spline and derivative, outlines 2d vectors left and right, 
# 1d vector angle assessments angle.

def semi_circle(N_angle):
    
    Radius=100
    Delta_radius=10
    
    
    #define track middle
    alfas= alfas = np.ones(30)*0.5
    theta = np.linspace(0, np.pi, int(3*len(alfas)/4))
    
    #inner and outer circles
    line1 = np.linspace(-Radius-Delta_radius+2*(Radius+Delta_radius)/
            (len(alfas)-len(theta)),Radius+Delta_radius-2*(Radius+Delta_radius)/
            (len(alfas)-len(theta)),int(len(alfas)-len(theta)))
    line2 = np.linspace(-Radius+Delta_radius+2*(Radius-Delta_radius)/
            (len(alfas)-len(theta)),Radius-Delta_radius-2*(Radius-Delta_radius)/
            (len(alfas)-len(theta)),int(len(alfas)-len(theta)))
    
    
    #outer line
    right = [np.concatenate(((Radius+Delta_radius)*np.cos(theta),line1)),
             np.concatenate(((Radius+Delta_radius)*np.sin(theta),
                             np.zeros(int(len(alfas)-len(theta)))))]
    
    #inner line
    left = [np.concatenate(((Radius-Delta_radius)*np.cos(theta),line2)),
            np.concatenate(((Radius-Delta_radius)*np.sin(theta),
                            np.zeros(int(len(alfas)-len(theta)))))]
    
    
    #calculate the spline (scipy spline), its derivative (scipy spline), 
    # angle (array over N_angle points)
    #left and right are the outlines and alfas is the random path
    spline, derivative, angle = path_info(left, right, alfas,N_angle)
    
    
    return spline, derivative, angle, right, left











#chooses the type of path. Calculates the scipy spline, scipy spline 
# derivative, assesses angles, define oulines left and right
#Inputs: name, external and internal .kml if existent, N_angles number 
#of assessment points for angles
#Output scipy spline and derivative, outlines 2d vectors left and right, 
# 1d vector angle assessments angle.

def choose_path(path_name,external,internal,N_angle):
    if path_name == "circle":
        spline, derivative, angle, right, left = circle(N_angle)
    elif path_name == "google_earth":      
        spline, derivative, angle, right, left = google_earth_path(external
                                                        ,internal,N_angle)
    elif path_name == "semi_circle":
        spline, derivative, angle, right, left = semi_circle(N_angle)
    else:
        print("Error: wrong path name")
        exit()
    return spline, derivative, angle, right, left
        
    

    
    