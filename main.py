##################################################################
###                            Import                          ###
##################################################################


#External libraries

import numpy as np
import faulthandler
faulthandler.enable()
import timeit
import matplotlib.pyplot as plt

#Internal functions

from Map_processing.choose_path import choose_path
from Physics.model1 import model1
from matplotlib.animation import FuncAnimation
from Physics.translate import translate_velocity, translate_acceleration
from Visualization.print import print_separator, print_table
from Simulation.optimization_main import *






##################################################################
###                     Problem definition                     ###
##################################################################


n_discretization=10 #number of path sections
N_path_points=1000 #plotting discretization
xsi = 1 #optimization scalar


#choose path
path_name = "circle"
external = 'Map_processing/Maps_kml/extHORTO.kml'
internal = 'Map_processing/Maps_kml/intHORTO.kml'



#create path splines and path spline derivative, assess orientation angles 
#over the sections, define outlines 
spline, derivative, angle, right, left = choose_path(path_name,external,
internal,N_angle=n_discretization)



#spline points for plotting
spline_points = spline(np.linspace(0,1,num = N_path_points))



#Define physics over the path
R_t, M_t, C_t, A_t = model1(spline,n_discretization)





##################################################################
###                           Choose Model                     ###
##################################################################


#Comment the models you dont want to compute

#Model abu
# t1_abu=init_optimization_abu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False) 

#Model bu
# t1_bu=init_optimization_bu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False) 

#Model b
# t1_b=init_optimization_b(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False)

#Model SOCP abu
# t1_SOCP_abu=init_optimization_SOCP_abu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False)

#Model SOCP b
t1_SOCP_b=init_optimization_SOCP_b(
    R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False)






##################################################################
###                 Model Performance Comparison              ###
##################################################################

#number of timeit assessments
N_computation_average=10


#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models = ["Time SOCP b"]


#Use same order as the models above
#t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]
results = [t1_SOCP_b[-1]]


#Call the timeit
model_performance(models,results,N_computation_average,R_t, M_t, C_t, 
     A_t,n_discretization,xsi,display=False)






##################################################################
###                     Real Path Calculation                  ###
##################################################################


#Only "Time abu" and "Time bu" available
controlled_path = controlled_path("Time bu",R_t, M_t, C_t, A_t,n_discretization,
                    xsi,spline_points,derivative,N_path_points)








##################################################################
###                     Circular Path test                     ###
##################################################################


#calculate absolute velocity
v = translate_velocity(derivative,decision_variables.x[0:n_discretization],
                       n_discretization)
#calculate absolute acceleration
a = translate_acceleration(derivative, derivative.derivative(),
                           decision_variables.x[0:n_discretization],
    decision_variables.x[n_discretization:2*n_discretization-1],
    n_discretization)

#Theoretical Circle maximum velocity
theoretical_v = np.ones(len(v))*np.sqrt(9.81*100)
#Theoretical Circle maximum acceleration
theoretical_a = 9.81*np.ones(len(v))

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8))
ax1.set_title('Velocity circular path')
ax1.set_xlabel('Path position')
ax1.set_ylabel('Velocity m/s')
ax2.set_title('Acceleration circular path')
ax2.set_xlabel('Path position')
ax2.set_ylabel('Acceleration m/sˆ2')
ax3.set_title('Work')
ax3.set_xlabel('Path position')
ax3.set_ylabel('Energy J')
ax1.set_ylim(0,max(v)+5)
ax2.set_ylim(0,max(a)+5)
ax3.set_ylim(min(E1)-20,max(E1)+20)
ax1.plot(np.linspace(0,1,n_discretization-1),
         v,'-r',label="Optimized velocity")
ax1.plot(np.linspace(0,1,n_discretization-1),
         theoretical_v,'-b',label="Theoretical velocity")
ax2.plot(np.linspace(0,1,n_discretization-1),
         a,'*r',label="Optimized acceleration")
ax2.plot(np.linspace(0,1,n_discretization-1),
         theoretical_a,'-b',label="Theoretical acceleration")
ax3.plot(np.linspace(0,1,n_discretization-1),E1,'-b',label="Work")
ax1.grid()
ax2.grid()
ax3.grid()
ax1.legend()
ax2.legend()
ax3.legend()
plt.show







##################################################################
###                        Animation Plot                      ###
##################################################################


#spline discretization over sections 
spline_points_animation = spline(np.linspace(0,1,num = n_discretization))


#create figure
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))


#Animation speed
animation_time = 100


# Set limits and labels for the first subplot
ax1.set_xlim(min(right[0])-10,max(right[0])+10)
ax1.set_ylim(min(right[1])-10,max(right[1])+10)
ax1.set_title('First guess constant dtheta Animation')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')


# Set limits and labels for the second subplot
ax2.set_xlim(min(right[0])-10,max(right[0])+10)
ax2.set_ylim(min(right[1])-10,max(right[1])+10)
ax2.set_title('Optimized Animation')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')


#Set the friction circle radius
radius = 85 * 9.81


# Create an array of angles for friction circle
theta = np.linspace(0, 2 * np.pi, 100)


# Parametric equations for the frictioncircle
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
#ax1.plot(controlled_path0[0], controlled_path0[1], 'k--', label='Real path')
#ax2.plot(controlled_path1[0], controlled_path1[1], 'k--', label='Real path')
ax2.plot(spline_points[0], spline_points[1], 'g--', label='Path')
ax2.plot(right[0], right[1], 'gray')
ax2.plot(left[0], left[1], 'red')
ax3.plot(x,y,'.r')
ax4.plot(x,y,'.r')

# Create lines for both objects
line1, = ax1.plot([], [], 'bo', label='Object 1')  # Blue point for object 1
line2, = ax2.plot([], [], 'ro', label='Object 2')  # Red point for object 2

# Create lines for the forces (bottom row)
force_line1, = ax3.plot([], [], 'b.', label='Force 1')  # Blue for obj 1 force
force_line2, = ax4.plot([], [], 'b.', label='Force 2')  # Red for obj 2 force

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
        line1.set_data([spline_points_animation[0, frame]], 
                       [spline_points_animation[1, frame]])
        force_line1.set_data(forcex0[:frame], forcey0[:frame]) #Update force data
    return line1, force_line1

def update2(frame):
    if frame < len(t1):
        line2.set_data([spline_points_animation[0, frame]], 
                    [spline_points_animation[1, frame]])
        force_line2.set_data(forcex1[:frame], forcey1[:frame]) #Update force data
    return line2, force_line2





def next_frame1(i):
    return update1(i)

def next_frame2(i):
    return update2(i)

# Create the animations
ani1 = FuncAnimation(fig, next_frame1, frames=len(t0), init_func=init1,
blit=True, repeat=True, interval=(t0[1] - t0[0]) * animation_time 
if len(t0) > 1 else animation_time)
ani2 = FuncAnimation(fig, next_frame2, frames=len(t1), init_func=init2,
blit=True, repeat=True, interval=(t1[1] - t1[0]) * animation_time 
if len(t1) > 1 else animation_time)

# Adjust the timing for each animation based on its time vector
 # Average time difference for animation 1
ani1.event_source.interval = np.mean(np.diff(t0)) * animation_time 
# Average time difference for animation 2
ani2.event_source.interval = np.mean(np.diff(t1)) * animation_time 

# Show the animation
plt.tight_layout()
plt.show()



##################################################################
###                    Comparison Methods                      ###
##################################################################


#Absolut velocity a,b,u method
v = translate_velocity(derivative,decision_variables.x[0:n_discretization],
                       n_discretization)

#Absolut velocity bu method
v_bu = translate_velocity(derivative,decision_variables_bu.x[0:n_discretization],
                          n_discretization)

#Absolut velocity b method
v_b = translate_velocity(derivative,decision_variables_b.x[0:n_discretization],
                         n_discretization)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
ax1.set_title('Comparison Velocities')
ax1.set_xlabel('Path position')
ax1.set_ylabel('Velocity m/s')
ax2.set_title('Time to traverse')
ax2.set_xlabel('Path position')
ax2.set_ylabel('Time')
ax1.set_ylim(0,max(max(v),max(v_b))+5)
ax2.set_ylim(0,max(max(t1),max(t1_bu))+5)
ax1.plot(np.linspace(0,1,n_discretization-1),v,'-r',
         label="Optimized velocity abu")
ax1.plot(np.linspace(0,1,n_discretization-1),v_bu,'*b',
         label="Optimized velocity bu")
ax1.plot(np.linspace(0,1,n_discretization-1),v_b,'og',
         label="Optimized velocity b")
ax2.plot(np.linspace(0,1,n_discretization),t1,'-r',label="Time abu")
ax2.plot(np.linspace(0,1,n_discretization),t1_bu,'*b',label="Time bu")
ax2.plot(np.linspace(0,1,n_discretization),t1_b,'og',label="Time b")
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()
plt.tight_layout()
plt.show()

