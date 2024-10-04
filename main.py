import numpy as np
#import pygame as gm
import matplotlib.pyplot as plt
from Map_processing.choose_path import choose_path
from Physics.model1 import model1
from Simulation.optimization import optimization
from Simulation.reconstruct import reconstruct, interpolate_u, control_system
from matplotlib.animation import FuncAnimation
from Physics.translate import translate_velocity, translate_acceleration


n_discretization=30 #number of path sections
N_path_points=1000 #plotting discretization

#choose path
path_name = "google_earth"
external = 'Map_processing/Maps_kml/extHORTO.kml'
internal = 'Map_processing/Maps_kml/intHORTO.kml'

#create path splines and path spline derivative, assess orientation angles over the sections, define outlines 
spline, derivative, angle, right, left = choose_path(path_name,external,internal,N_angle=n_discretization)
spline_points = spline(np.linspace(0,1,num = N_path_points))

#Define physics over the path
R_t, M_t, C_t, A_t = model1(spline,n_discretization)


xsi = 1 #optimization scalar

#finds the optimal solution and innitial guess. Outputs generalized velocity square b, generalized acceleration a, and forces u
decision_variables, x0 = optimization(R_t, M_t, C_t, A_t,n_discretization,xsi)

#extract the forces from the flattened result array
forcex0=x0[2*n_discretization-1:3*n_discretization-2]
forcey0=x0[3*n_discretization-2:len(decision_variables.x)]
forcex1=decision_variables.x[2*n_discretization-1:3*n_discretization-2]
forcey1=decision_variables.x[3*n_discretization-2:len(decision_variables.x)]

#calculated time to run each trajectory using generalized velocity square b 
t0 = reconstruct(x0[0:n_discretization])
t1=reconstruct(decision_variables.x[0:n_discretization])

#calculate command vector if needed
command_vector = interpolate_u(np.transpose([forcex1,forcey1]),t1,num_points=1000)

#calculate initial position and velocity for innitial guess 
x00 =[spline_points[0][0], spline_points[1][0]]
v00 = [derivative([0])[0][0]*np.sqrt(x0[0]+x0[1]/2),derivative([0])[1][0]*np.sqrt(x0[0]+x0[1]/2)]

#Calculates the real path with the control of initial guess
controlled_path0 = control_system([forcex0, forcey0],x00,v00,t0,N_path_points)

#calculate initial position and velocity for optimized trajectory 
x10 =[spline_points[0][0], spline_points[1][0]]
v10 = [derivative([0])[0][0]*np.sqrt(decision_variables.x[0]+decision_variables.x[1]/2),derivative([0])[1][0]*np.sqrt(decision_variables.x[0]+decision_variables.x[1]/2)]
#Calculates the real path with the control of optimized trajectory
controlled_path1 = control_system([forcex1, forcey1],x10,v10,t1,N_path_points)


#calculate absolute velocity
v = translate_velocity(derivative,decision_variables.x[0:n_discretization],n_discretization)
#calculate absolute acceleration
a = translate_acceleration(derivative, derivative.derivative(),decision_variables.x[0:n_discretization], decision_variables.x[n_discretization:2*n_discretization-1],n_discretization)

#Theoretical Circle maximum velocity
theoretical_v = np.ones(len(v))*np.sqrt(9.81*100)
#Theoretical Circle maximum acceleration
theoretical_a = 9.81*np.ones(len(v))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
ax1.set_title('Velocity circular path')
ax1.set_xlabel('Path position')
ax1.set_ylabel('Velocity m/s')
ax2.set_title('Acceleration circular path')
ax2.set_xlabel('Path position')
ax2.set_ylabel('Acceleration m/sˆ2')
ax1.set_ylim(0,max(v)+5)
ax2.set_ylim(0,max(a)+5)
ax1.plot(np.linspace(0,1,n_discretization-1),v,'*r',label="Optimized velocity")
ax1.plot(np.linspace(0,1,n_discretization-1),theoretical_v,'-b',label="Theoretical velocity")
ax2.plot(np.linspace(0,1,n_discretization-1),a,'*r',label="Optimized acceleration")
ax2.plot(np.linspace(0,1,n_discretization-1),theoretical_a,'-b',label="Theoretical acceleration")
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()
plt.show

# fig = plt.figure()
# plt.plot(spline_points[0],spline_points[1])
# plt.plot(right[0], right[1], 'gray')
# plt.plot(left[0], left[1], 'red')
# plt.show

# fig = plt.figure()
# plt.plot(np.linspace(0,1,num = N),angle*360/2/3.14159) 
# plt.plot(np.linspace(0,1,num = N),np.sin(angle)*100)
# plt.show()



#animation

spline_points_animation = spline(np.linspace(0,1,num = n_discretization))
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))


animation_time = 100
# Set limits and labels for the first subplot
ax1.set_xlim(min(min(controlled_path0[0]),min(right[0]))-10,max(max(controlled_path0[0]),max(right[0]))+10)
ax1.set_ylim(min(min(controlled_path0[1]),min(right[1]))-10,max(max(controlled_path0[1]),max(right[1]))+10)
ax1.set_title('First guess constant dtheta Animation')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

# Set limits and labels for the second subplot
ax2.set_xlim(min(min(controlled_path1[0]),min(right[0]))-10,max(max(controlled_path1[0]),max(right[0]))+10)
ax2.set_ylim(min(min(controlled_path1[1]),min(right[1]))-10,max(max(controlled_path1[1]),max(right[1]))+10)
ax2.set_title('Optimized Animation')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')

radius = 85 * 9.81

# Create an array of angles
theta = np.linspace(0, 2 * np.pi, 100)

# Parametric equations for the circle
x = radius * np.cos(theta)
y = radius * np.sin(theta)


# Set limits and labels for the force plot 1 (bottom-left)
ax3.set_xlim(-85*9.81, 85*9.81)  # X-axis is time
ax3.set_ylim(-85*9.81, 85*9.81)    # Y-axis is force
ax3.set_title('First guess Force')
ax3.set_xlabel('Time')
ax3.set_ylabel('Force')
ax3.grid()

# Set limits and labels for the force plot 2 (bottom-right)
ax4.set_xlim(-85*9.81, 85*9.81)    # X-axis is time
ax4.set_ylim(-85*9.81, 85*9.81)    # Y-axis is force
ax4.set_title('Optimized Force')
ax4.set_xlabel('Time')
ax4.set_ylabel('Force')
ax4.grid()




# Plot the spline points in both axes
ax1.plot(spline_points[0], spline_points[1], 'g--', label='Path')
ax1.plot(right[0], right[1], 'gray')
ax1.plot(left[0], left[1], 'red')
ax1.plot(controlled_path0[0], controlled_path0[1], 'k--', label='Real path')
ax2.plot(controlled_path1[0], controlled_path1[1], 'k--', label='Real path')
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
ani1 = FuncAnimation(fig, next_frame1, frames=len(t0), init_func=init1, blit=True, repeat=True, interval=(t0[1] - t0[0]) * animation_time if len(t0) > 1 else animation_time)
ani2 = FuncAnimation(fig, next_frame2, frames=len(t1), init_func=init2, blit=True, repeat=True, interval=(t1[1] - t1[0]) * animation_time if len(t1) > 1 else animation_time)

# Adjust the timing for each animation based on its time vector
ani1.event_source.interval = np.mean(np.diff(t0)) * animation_time  # Average time difference for animation 1
ani2.event_source.interval = np.mean(np.diff(t1)) * animation_time # Average time difference for animation 2

# Show the animation
plt.tight_layout()
plt.show()







