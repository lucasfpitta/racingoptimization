import numpy as np
#import pygame as gm
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import pandas as pd
import scipy as scp
import random as rnd
from Map_processing.mapreader import get_outline
from Map_processing.point_reduction_and_smooting import reduce_points, smooth_moving_average
from splines.splines import path_info
from Physics.model1 import model1
from Simulation.optimization import optimization
from Simulation.reconstruct import reconstruct
from matplotlib.animation import FuncAnimation

alfa = 0.01
n_neighborhood = 16
max_slope = 0.06
N=1000

#read the map and set the metric coordinates
right, left,  = get_outline('Map_processing/Maps_kml/extHORTO.kml','Map_processing/Maps_kml/intHORTO.kml')
coords= reduce_points(right,left)
#smooth array 2 and 5 (z coordinate)
#coords[2], coords[5] = smooth_moving_average(coords[0:3], coords[3:6],int(len(coords[0])/n_neighborhood),alfa,max_slope)


"""fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(right[0], right[1], right[2], 'gray')
ax.plot3D(left[0], left[1], left[2], 'red')
plt.show()"""

alfas = np.random.random_sample(len(right[0]))
alfas[-1]=alfas[0]
spline, derivative, angle = path_info(left, right, alfas,N)
spline_points = spline(np.linspace(0,1,num = N))
n_discretization=50
R_t, M_t, C_t, A_t = model1(spline,n_discretization)


xsi = 1

decision_variables, x0 = optimization(R_t, M_t, C_t, A_t,n_discretization,xsi)
forcex0=x0[2*n_discretization-1:3*n_discretization-2]
forcey0=x0[3*n_discretization-2:len(decision_variables.x)]
forcex1=decision_variables.x[2*n_discretization-1:3*n_discretization-2]
forcey1=decision_variables.x[3*n_discretization-2:len(decision_variables.x)]

t0 = reconstruct(x0[0:n_discretization])
t1=reconstruct(decision_variables.x[0:n_discretization])

print(t0,t1)

# fig = plt.figure()
# plt.plot(spline_points[0],spline_points[1])
# plt.plot(right[0], right[1], 'gray')
# plt.plot(left[0], left[1], 'red')
# plt.show

# fig = plt.figure()
# plt.plot(np.linspace(0,1,num = N),angle*360/2/3.14159) 
# plt.plot(np.linspace(0,1,num = N),np.sin(angle)*100)
# plt.show()





spline_points_animation = spline(np.linspace(0,1,num = n_discretization))
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))

# Set limits and labels for the first subplot
ax1.set_xlim(-100, 70)
ax1.set_ylim(-50, 600)
ax1.set_title('First guess constant dtheta Animation')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

# Set limits and labels for the second subplot
ax2.set_xlim(-100, 70)
ax2.set_ylim(-50, 600)
ax2.set_title('Optimized Animation')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')

radius = 30 * 9.81

# Create an array of angles
theta = np.linspace(0, 2 * np.pi, 100)

# Parametric equations for the circle
x = radius * np.cos(theta)
y = radius * np.sin(theta)


# Set limits and labels for the force plot 1 (bottom-left)
ax3.set_xlim(-30*9.81, 30*9.81)  # X-axis is time
ax3.set_ylim(-30*9.81, 30*9.81)    # Y-axis is force
ax3.set_title('First guess 1 Force')
ax3.set_xlabel('Time')
ax3.set_ylabel('Force')

# Set limits and labels for the force plot 2 (bottom-right)
ax4.set_xlim(-30*9.81, 30*9.81)    # X-axis is time
ax4.set_ylim(-30*9.81, 30*9.81)    # Y-axis is force
ax4.set_title('Optimized Force')
ax4.set_xlabel('Time')
ax4.set_ylabel('Force')




# Plot the spline points in both axes
ax1.plot(spline_points[0], spline_points[1], 'g--', label='Path')
ax1.plot(right[0], right[1], 'gray')
ax1.plot(left[0], left[1], 'red')
ax2.plot(spline_points[0], spline_points[1], 'g--', label='Path')
ax2.plot(right[0], right[1], 'gray')
ax2.plot(left[0], left[1], 'red')
ax3.plot(x,y,'.r')
ax4.plot(x,y,'.r')

# Create lines for both objects
line1, = ax1.plot([], [], 'bo', label='Object 1')  # Blue point for object 1
line2, = ax2.plot([], [], 'ro', label='Object 2')  # Red point for object 2

# Create lines for the forces (bottom row)
force_line1, = ax3.plot([], [], 'b.', label='Force 1')  # Blue line for object 1 force
force_line2, = ax4.plot([], [], 'b.', label='Force 2')  # Red line for object 2 force

# Add legends to all axes
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()


def init1():
    line1.set_data([], [])
    force_line1.set_data([], [])
    return line1, force_line1

def init2():
    line2.set_data([], [])
    force_line2.set_data([], [])
    return line2, force_line2

# Update function for object animation
def update1(frame):
    if frame < len(t0):
        line1.set_data([spline_points_animation[0, frame]], [spline_points_animation[1, frame]])
        force_line1.set_data(forcex0[:frame], forcey0[:frame])  # Update force data
    return line1, force_line1

def update2(frame):
    if frame < len(t1):
        line2.set_data([spline_points_animation[0, frame]], [spline_points_animation[1, frame]])
        force_line2.set_data(forcex1[:frame], forcey1[:frame])  # Update force data
    return line2, force_line2





def next_frame1(i):
    return update1(i)

def next_frame2(i):
    return update2(i)

# Create the animations
ani1 = FuncAnimation(fig, next_frame1, frames=len(t0), init_func=init1, blit=True, repeat=True, interval=(t0[1] - t0[0]) * 1000 if len(t0) > 1 else 1000)
ani2 = FuncAnimation(fig, next_frame2, frames=len(t1), init_func=init2, blit=True, repeat=True, interval=(t1[1] - t1[0]) * 1000 if len(t1) > 1 else 1000)

# Adjust the timing for each animation based on its time vector
ani1.event_source.interval = np.mean(np.diff(t0)) * 1000  # Average time difference for animation 1
ani2.event_source.interval = np.mean(np.diff(t1)) * 1000  # Average time difference for animation 2

# Show the animation
plt.tight_layout()
plt.show()







