from Map_processing.mapreader import get_outline
from Map_processing.point_reduction_and_smooting import reduce_points, smooth_moving_average
from splines.splines import path_info
import numpy as np


def google_earth_path(external,internal,N_angle):
    alfa = 0.01
    n_neighborhood = 16
    max_slope = 0.06
    #read the map and set the metric coordinates
    right, left,  = get_outline(external,internal)#('Map_processing/Maps_kml/extHORTO.kml','Map_processing/Maps_kml/intHORTO.kml')
    #coords= reduce_points(right,left)
    #smooth array 2 and 5 (z coordinate)
    #coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],int(len(coords[0])/n_neighborhood),alfa,max_slope)
    
    alfas = np.random.random_sample(len(right[0]))
    alfas[-1]=alfas[0]
    
    spline, derivative, angle = path_info(left, right, alfas,N_angle)
    return spline, derivative, angle, right, left

def circle(N_angle):
    Radius=100
    Delta_radius=10
    alfas= alfas = np.ones(30)*0.5
    theta = np.linspace(0, 2 * np.pi, len(alfas))
    right = [(Radius+Delta_radius)*np.cos(theta),(Radius+Delta_radius)*np.sin(theta)]
    left = [(Radius-Delta_radius)*np.cos(theta),(Radius-Delta_radius)*np.sin(theta)]
    
    spline, derivative, angle = path_info(left, right, alfas,N_angle)
    return spline, derivative, angle, right, left

# def semi_circle(ext,int,alfas,N):
    
    
#     return spline, derivative, angle, right, left

def choose_path(path_name,external,internal,N_angle):
    if path_name == "circle":
        spline, derivative, angle, right, left = circle(N_angle)
    elif path_name == "google_earth":      
        spline, derivative, angle, right, left = google_earth_path(external,internal,N_angle)
    else:
        print("Error: wrong path name")
        exit()
    return spline, derivative, angle, right, left
        
    

    
    