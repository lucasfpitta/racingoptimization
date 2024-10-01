import numpy as np
import pygame as gm
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import pandas as pd
import scipy as scp
import random as rnd
from Map_processing.mapreader import get_outline
from Map_processing.point_reduction_and_smooting import reduce_points, smooth_moving_average
from splines.splines import path_info
from Physics.model1 import model1

alfa = 0.01
n_neighborhood = 16
max_slope = 0.06
N=1000

#read the map and set the metric coordinates
right, left,  = get_outline('Map_processing/Maps_kml/extHORTO.kml','Map_processing/Maps_kml/intHORTO.kml')
coords= reduce_points(right,left)
#smooth array 2 and 5 (z coordinate)
#coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],int(len(coords[0])/n_neighborhood),alfa,max_slope)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(right[0], right[1], right[2], 'gray')
ax.plot3D(left[0], left[1], left[2], 'red')
plt.show()

alfas = np.random.random_sample(len(right[0]))
alfas[-1]=alfas[0]
spline, derivative, angle = path_info(left, right, alfas,N)
spline_points = spline(np.linspace(0,1,num = N))
M=10
R_t, M_t, C_t, A_t = model1(spline,M)
print("A_t",A_t)




fig = plt.figure()
plt.plot(spline_points[0],spline_points[1])
plt.plot(right[0], right[1], 'gray')
plt.plot(left[0], left[1], 'red')
plt.show

fig = plt.figure()
plt.plot(np.linspace(0,1,num = N),angle*360/2/3.14159) 
plt.plot(np.linspace(0,1,num = N),np.sin(angle)*100)
plt.show()


