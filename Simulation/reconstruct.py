import numpy as np


#Calculates the time to traverse each section using the generalized velocity squared b
#Input 1d b vector
#Output 1d time vector
def reconstruct(b):
    dimension = len(b)
    t=np.zeros(dimension)
    for i in range(dimension-1):
        t[i+1] = 2*1/(dimension-1)/(b[i]**0.5+b[i+1]**0.5)+t[i]
    return t

#Calulates U vector over a evenly space time vector for future control
#Input 1d b vector
#Output 1d time vector

def interpolate_u(u, t1, num_points):
    # Create the new evenly spaced time vector
    new_t = np.linspace(0, t1[-1], num_points)
    
    # Initialize new u vector (same length as new_t)
    new_u = np.zeros((num_points,2))  
    
    # Loop through the new_t vector and assign values based on the original t1 intervals
    for i in range(len(t1) - 1):
        # Find the indices in new_t that fall within the current interval [t1[i], t1[i+1]]
        index = (new_t >= t1[i]) & (new_t < t1[i+1])
        
        # Assign the corresponding u value to new_u for this interval
        new_u[index] = u[i]

    # Handle the last interval edge case
    new_u[new_t >= t1[-1]] = u[-1]

    return new_u

    # Approximates the force for a given time t for the real path calculation based on the control
    #Input the time t used by the integrator, the vector of time to traverse each section t1, the corresponding 2d vector u for each section
    # Output 2d vector of force at time t.

def approximate_force(t, t1, u):
    # Find the first index where t1[i] is greater than t
    idx = np.searchsorted(t1, t, side='left')-1
    if idx < len(t1):
        # Return the corresponding force value given that the force is constant in the interval
        return [u[0][idx], u[1][idx]]
    else:
        # If t is beyond the last time point, return the last force value 
        return u[-1]
    
import numpy as np
from scipy.integrate import solve_ivp


#ODE solver to calculate real path
#Inputs 2d vector of forces in each section, initial position and initial velocity, time to traverse each section, discretization for the ODE
def control_system(u1,x0,v0,t1,n_discretization):
    #Ode function
    #Input from the solver scalar time t and state vector z
    #Other inputs as arguments, matrices R, M, C from the physical model, force in each section and time to traverse each section
    #Output derivatives
    def ode_system(t, z, R, M, C, u1, t1):
        x = z[:len(z)//2]  # First half of the state vector is x
        v = z[len(z)//2:]  # Second half is v 

        # First-order system
        dxdt = v
        #to calculate dv/dt we find u(t)
        u=approximate_force(t,t1,u1)
        #dv/dt comes from the state equation
        dvdt = np.linalg.inv(M) @ (R @ u - C @ v)

        return np.concatenate([dxdt, dvdt])

    # Physics setup
    R = np.identity(2)  #R matrix
    M = np.identity(2)*85  #M matrix
    C = np.identity(2)*0 #C matrix

    #initial state vector
    z0 = np.concatenate([x0, v0])

    #Time span for the integration
    t_span = (0, t1[-1])
    t_eval = np.linspace(0, t1[-1], n_discretization)

    #Solve the system
    sol = solve_ivp(ode_system, t_span, z0, args=(R, M, C, u1, t1), t_eval=t_eval)

    # Extract solution for x and v
    x_sol = sol.y[:len(x0), :]
    return x_sol



