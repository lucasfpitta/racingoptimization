##################################################################
###                            Import                          ###
##################################################################


#External libraries

import numpy as np
import timeit



#Internal functions
from Simulation.optimization_abu import optimization_abu
from Simulation.optimization_b import optimization_b
from Simulation.optimization_bu import optimization_bu
from Simulation.optimization_SOCP_abu import optimization_SOCP_abu
from Simulation.optimization_SOCP_b import optimization_SOCP_b
from Simulation.reconstruct import reconstruct, interpolate_u, control_system
from Visualization.print import print_separator, print_table








##################################################################
### Velocity, acceleration and force optimization Model (abu)  ###
##################################################################



def init_optimization_abu(R_t, M_t, C_t, A_t,n_discretization,
                                      xsi,display,plot):
    
     if display: 
          print_separator(
               "Velocity, acceleration and force optimization Model (abu)")

     #finds the optimal solution and innitial guess. Outputs generalized 
     #velocity square b, 
     #generalized acceleration a, and forces u
     decision_variables_abu, x0_abu = optimization_abu(R_t, M_t, C_t, A_t,
                                             n_discretization,xsi,display)


     #extract the forces from the flattened result array
     forcex0_abu=x0_abu[2*n_discretization-1:3*n_discretization-2]
     forcey0_abu=x0_abu[3*n_discretization-2:len(decision_variables_abu.x)]
     forcex1_abu=decision_variables_abu.x[2*n_discretization-1:3*
                                                       n_discretization-2]
     forcey1_abu=decision_variables_abu.x[3*n_discretization-2:
                                             len(decision_variables_abu.x)]



     #calculated time to run each trajectory using generalized velocity square b 
     t0_abu = reconstruct(x0_abu[0:n_discretization])
     t1_abu=reconstruct(decision_variables_abu.x[0:n_discretization])



     if plot:
          return t0_abu,t1_abu,forcex0_abu,forcey0_abu,forcex1_abu,\
               forcey1_abu,x0_abu,decision_variables_abu.x
     else:
          return t1_abu







##################################################################
###              Velocity Force optimization Model (bu)        ###
##################################################################


def init_optimization_bu(R_t, M_t, C_t, A_t,n_discretization,
                                      xsi,display,plot):
    
     if display:
          print_separator("Velocity Force optimization Model (bu)")

     #finds the optimal solution and innitial guess. Outputs generalized 
     #velocity square b
     decision_variables_bu, x0_bu = optimization_bu(R_t, M_t, C_t, A_t,
                                             n_discretization,xsi,display)


     #calculated time to run each trajectory using generalized velocity 
     #square b 
     t0_bu = reconstruct(x0_bu[0:n_discretization])
     t1_bu=reconstruct(decision_variables_bu.x[0:n_discretization])


     #extract the forces from the flattened result array
     forcex0_bu=x0_bu[n_discretization:2*n_discretization-1]
     forcey0_bu=x0_bu[2*n_discretization-1:len(decision_variables_bu.x)]
     forcex1_bu=decision_variables_bu.x[n_discretization:2*n_discretization-1]
     forcey1_bu=decision_variables_bu.x[2*n_discretization-1:
          len(decision_variables_bu.x)]
     
     if plot:
          return t0_bu,t1_bu,forcex0_bu,forcey0_bu,forcex1_bu,\
               forcey1_bu,x0_bu,decision_variables_bu.x
     else:
          return t1_bu







##################################################################
###                Velocity optimization Model (b)             ###
##################################################################


def init_optimization_b(R_t, M_t, C_t, A_t,n_discretization,
                                      xsi,display,plot):

     if display:
          print_separator("Velocity optimization Model (b)")

     #finds the optimal solution and innitial guess. Outputs generalized 
     #velocity square b
     decision_variables_b, x0_b = optimization_b(R_t, M_t, C_t, A_t,
                                                  n_discretization,xsi,display)


     #calculated time to run each trajectory using generalized velocity 
     #square b 
     t0_b = reconstruct(x0_b[0:n_discretization])
     t1_b=reconstruct(decision_variables_b.x[0:n_discretization])
     
     if plot:
          return t0_b,t1_b,x0_b,decision_variables_b.x
     else:
          return t1_b








##################################################################
###                 Second-order Cone (abu) Model              ###
##################################################################

def init_optimization_SOCP_abu(R_t, M_t, C_t, A_t,n_discretization,
                                      xsi,display,plot):
     if display:
          print_separator("Second-order Cone (abu) Model")

     #finds the optimal solution. Outputs vector with variables a, b, u, 
     # c, d 
     decision_variables_SOCP_abu = optimization_SOCP_abu(R_t, M_t, C_t, 
                                             A_t,n_discretization,xsi,display)

     #calculated time to run each trajectory using generalized velocity 
     #square b 
     t1_SOCP_abu=reconstruct(decision_variables_SOCP_abu[n_discretization-1
                                                       :2*n_discretization-1])
     
     if plot:
          return t1_SOCP_abu,decision_variables_SOCP_abu
     else:
          return t1_SOCP_abu







##################################################################
###                  Second-order Cone (b) Model               ###
##################################################################

def init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,
                                      xsi,display,plot):
     if display:
          print_separator("Second-order Cone (b) Model")

     #finds the optimal solution. Outputs vector with variables a, b, u, 
     # c, d 
     decision_variables_SOCP_b = optimization_SOCP_b(R_t, M_t, C_t, 
                                             A_t,n_discretization,xsi,display)


     #calculated time to run each trajectory using generalized velocity 
     #square b 
     t1_SOCP_b=reconstruct(decision_variables_SOCP_b[0:n_discretization])
     
     if plot:
          return t1_SOCP_b,decision_variables_SOCP_b
     else:
          return t1_SOCP_b
   
   
   
   
   


##################################################################
###                 Model Performance Comparison              ###
##################################################################

def model_performance(models,results,N_computation_average,R_t, M_t, C_t, 
     A_t,n_discretization,xsi,display):
     print_separator("Model Performance Comparison")
     print(f"Number of sections: {n_discretization}")
     
     
     if len(models)!=len(results):
          print("There's a diferent number of models and results")
          SystemExit
     if len(models)==0:
          print("No time assessment requested")
          return
     
     
     # Dictionary of Models
     Models_dict = {
          "Time abu": init_optimization_abu,
          "Time bu": init_optimization_bu,
          "Time b": init_optimization_b,
          "Time SOCP abu": init_optimization_SOCP_abu,
          "Time SOCP b": init_optimization_SOCP_b
          }
     
     #Dictionary with the times to compute
     compute_times_dict = {}
     
     computation_time = []
     for name, func in Models_dict.items():
          time_taken = timeit.timeit(lambda:func(R_t, M_t, C_t, 
          A_t,n_discretization,xsi,display,plot=False), 
                              number=N_computation_average)
          compute_times_dict[name] = time_taken


     computation_time = [compute_times_dict[name] for name in models]
     print_table(models,results,computation_time)
     return computation_time







##################################################################
###                     Real Path Calculation                  ###
##################################################################


def controlled_path(model,R_t, M_t, C_t, A_t,n_discretization,
                    xsi,spline_points,derivative,N_path_points):
     
     print_separator("Real Path Calculation")
     
     
     Models_dict = {
          "Time abu": init_optimization_abu,
          "Time bu": init_optimization_bu,
          #"Time b": init_optimization_b,
          #"Time SOCP abu": init_optimization_SOCP_abu,
          #"Time SOCP b": init_optimization_SOCP_b
               }
     
     if model in Models_dict:
          t0,t1,forcex0,forcey0,forcex1,forcey1,x0,decision_variables=\
               Models_dict[model](R_t, M_t, C_t, A_t,n_discretization,xsi,\
                    display=False,plot=True)
     else:
          print("Invalid model name")
          return

     #calculate command vector if needed
     command_vector = interpolate_u(np.transpose([forcex1,forcey1]),t1,
                                   num_points=1000)


     #calculate initial position and velocity for innitial guess 
     x00 =[spline_points[0][0], spline_points[1][0]]
     v00 = [derivative([0])[0][0]*np.sqrt(x0[0]+x0[1]/2),derivative([0])[1][0]
          *np.sqrt(x0[0]+x0[1]/2)]


     #Calculates the real path with the control of initial guess
     controlled_path0 = control_system([forcex0, forcey0],x00,v00,
                                       t0,N_path_points)


     #calculate initial position and velocity for optimized trajectory 
     x10 =[spline_points[0][0], spline_points[1][0]]
     v10 = [derivative([0])[0][0]*np.sqrt(decision_variables[0]+
     decision_variables[1]/2),derivative([0])[1][0]*np.sqrt(
          decision_variables[0]+decision_variables[1]/2)]


     #Calculates the real path with the control of optimized trajectory
     controlled_path1 = control_system([forcex1, forcey1],x10,v10,t1,
                                       N_path_points)


     print("Real Path calculated successfully")
     return controlled_path
