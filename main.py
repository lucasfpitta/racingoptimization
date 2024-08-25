import numpy as np
import pygame as gm
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import pandas as pd
import scipy as sci
from Map_processing.mapreader import get_outline
from Map_processing.point_reduction_and_smooting import reduce_points, smooth_moving_average

alfa = 0.01
n_neighborhood = 16
max_slope = 0.06

#read the map and set the metric coordinates
right, left,  = get_outline('Map_processing/Maps_kml/extSEM.kml','Map_processing/Maps_kml/intSEM.kml')
coords= reduce_points(right,left)
#smooth array 2 and 5 (z coordinate)
coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],int(len(coords[0])/n_neighborhood),alfa,max_slope)
plt.plot(right[0,:100], right[1,:100], 'gray',left[0,:100], left[1,:100], 'red')
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(right[0], right[1], right[2], 'gray')
ax.plot3D(left[0], left[1], left[2], 'red')
plt.show()
print("Hello World")