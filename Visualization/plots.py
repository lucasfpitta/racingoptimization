from Physics.translate import translate_velocity, translate_acceleration
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from Simulation.optimization_main import *






##################################################################
###                     Circular Path test                     ###
##################################################################

def circular_path_test(derivative,decision_variables,n_discretization):
    
    #calculate absolute velocity
    v = translate_velocity(derivative,decision_variables[0:n_discretization],
                        n_discretization)
    
    #calculate absolute acceleration
    a = translate_acceleration(derivative, derivative.derivative(),
                            decision_variables[0:n_discretization],
        decision_variables[n_discretization:2*n_discretization-1],
        n_discretization)


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
    ax1.plot(np.linspace(0,1,n_discretization-1),
            v,'-r',label="Optimized velocity")
    ax1.plot(np.linspace(0,1,n_discretization-1),
            theoretical_v,'-b',label="Theoretical velocity")
    ax2.plot(np.linspace(0,1,n_discretization-1),
            a,'*r',label="Optimized acceleration")
    ax2.plot(np.linspace(0,1,n_discretization-1),
            theoretical_a,'-b',label="Theoretical acceleration")
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    plt.show()
    











##################################################################
###                        Animation Plot                      ###
##################################################################


def animation_(spline,right,left,spline_points,forcex0,forcey0,forcex1,forcey1
               ,t0,t1,n_discretization):
    
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

def comparison_plot(derivative,R_t, M_t, C_t, A_t,n_discretization,xsi):

    #all results needed
    t0_abu,t1_abu,forcex0_abu,forcey0_abu,forcex1_abu,forcey1_abu,x0_abu,\
        decision_variables_abu = init_optimization_abu(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)
    
    t0_bu,t1_bu,forcex0_bu,forcey0_bu,forcex1_bu,forcey1_bu,x0_bu,\
        decision_variables_bu = init_optimization_bu(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)
    
    t0_b,t1_b,x0_b,decision_variables_b = init_optimization_b(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)
    
    t1_SOCP_abu,decision_variables_SOCP_abu = init_optimization_SOCP_abu(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)
    
    t1_SOCP_b,decision_variables_SOCP_b = init_optimization_SOCP_b(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)



    #Absolut velocity a,b,u method
    v = translate_velocity(derivative,decision_variables_abu[0:n_discretization],
                        n_discretization)

    #Absolut velocity bu method
    v_bu = translate_velocity(derivative,decision_variables_bu[0:n_discretization],
                            n_discretization)

    #Absolut velocity b method
    v_b = translate_velocity(derivative,decision_variables_b[0:n_discretization],
                            n_discretization)
    
    #Absolut velocity SOCP_abu method
    v_SOCP_abu = translate_velocity(derivative,decision_variables_SOCP_abu[\
        n_discretization-1:2*n_discretization-1],
                            n_discretization)
    
    #Absolut velocity b method
    v_SOCP_b = translate_velocity(derivative,decision_variables_SOCP_b[0:n_discretization],
                            n_discretization)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    ax1.set_title('Comparison Velocities')
    ax1.set_xlabel('Path position')
    ax1.set_ylabel('Velocity m/s')
    ax2.set_title('Time to traverse')
    ax2.set_xlabel('Path position')
    ax2.set_ylabel('Time')
    ax1.set_ylim(0,max(max(v),max(v_b))+5)
    ax2.set_ylim(0,max(max(t1_abu),max(t1_bu))+5)
    ax1.plot(np.linspace(0,1,n_discretization-1),v,'-r',
            label="Optimized velocity abu")
    ax1.plot(np.linspace(0,1,n_discretization-1),v_bu,'*b',
            label="Optimized velocity bu")
    ax1.plot(np.linspace(0,1,n_discretization-1),v_b,'og',
            label="Optimized velocity b")
    ax1.plot(np.linspace(0,1,n_discretization-1),v_SOCP_abu,'oy',
            label="Optimized velocity SOCP abu")
    ax1.plot(np.linspace(0,1,n_discretization-1),v_SOCP_b,'-k',
            label="Optimized velocity SOCP b")
    ax2.plot(np.linspace(0,1,n_discretization),t1_abu,'-r',label="Time abu")
    ax2.plot(np.linspace(0,1,n_discretization),t1_bu,'*b',label="Time bu")
    ax2.plot(np.linspace(0,1,n_discretization),t1_b,'og',label="Time b")
    ax2.plot(np.linspace(0,1,n_discretization),t1_SOCP_abu,'oy',label="Time SOCP abu")
    ax2.plot(np.linspace(0,1,n_discretization),t1_SOCP_b,'-k',label="Time SOCP b")
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    plt.tight_layout()
    plt.show()