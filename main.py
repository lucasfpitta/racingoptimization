##################################################################
###                            Import                          ###
##################################################################


#External libraries

import numpy as np
import faulthandler
faulthandler.enable()

#Internal functions

from Map_processing.choose_path import choose_path
from Physics.model1 import model1
from Simulation.optimization_main import *
from Visualization.plots import *
from Comparison.Opt_models_comparison import *







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
t1_abu=init_optimization_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False) 

#Model bu
t1_bu=init_optimization_bu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False) 

#Model b
t1_b=init_optimization_b(
    R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False)

#Model SOCP abu
t1_SOCP_abu=init_optimization_SOCP_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,display=True,plot=False)

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

models = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]


#Use same order as the models above
#t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]
results = [t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]]


#Call the timeit
computation_time = model_performance(models,results,N_computation_average,
                    R_t, M_t, C_t,A_t,n_discretization,xsi,display=False)



##################################################################
###                       Comparison Export                   ###
##################################################################

#number of sections to access
discretizations = [10,18,32,56,100]

#number of timeit assessments
N_computation_average=200

#chose the filename
filename = "Comparison/Results/comparison_timeit.csv"

#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models_export = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]

export_comparison_to_csv(models_export, discretizations,filename,
                         N_computation_average,xsi,spline)





#data Dictionary
data = read_csv_to_dict(filename)

#Complexity calculation

complexity = fit_log(data)

model_complexity(models,complexity)


##################################################################
###                     Real Path Calculation                  ###
##################################################################


#Only "Time abu" and "Time bu" available
controlled_path = controlled_path("Time bu",R_t, M_t, C_t, A_t,
        n_discretization,xsi,spline_points,derivative,N_path_points)









##################################################################
###                            Plots                           ###
##################################################################


#Uncomment the plots you want

#solution general model
t0_abu,t1_abu,forcex0_abu,forcey0_abu,forcex1_abu,forcey1_abu,x0_abu,\
    decision_variables_abu = init_optimization_abu(
          R_t, M_t, C_t, A_t,n_discretization,xsi,display=False,plot=True)


#Test if the circular path velocity is equal to the theoretical
circular_path_test(derivative,decision_variables_abu,n_discretization)

#Animates initial guess vs optimized solution
animation_(spline,right,left,spline_points,forcex0_abu,forcey0_abu,\
            forcex1_abu,forcey1_abu
               ,t0_abu,t1_abu,n_discretization)

#Solution comparison plot
comparison_plot(derivative,R_t, M_t, C_t, A_t,n_discretization,xsi)