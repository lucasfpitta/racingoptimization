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
right, left,  = get_outline('Map_processing/Maps_kml/catalunya_2022_right.kml','Map_processing/Maps_kml/catalunya_2022_left.kml')
coords= reduce_points(right,left)
#smooth array 2 and 5 (z coordinate)
coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],int(len(coords[0])/n_neighborhood),alfa,max_slope)