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
from Physics.model2 import model2
from Physics.model3 import model3
from Physics.model4 import model4
from Simulation.optimization_main import *
from Visualization.plots import *
from Comparison.Opt_models_comparison import *
from splines.splines import model4_extra_angles





##################################################################
###                     Problem definition                     ###
##################################################################


#Optimization variables
n_discretization=15 #number of path sections
N_path_points=1000 #plotting discretization
xsi = 1 #optimization scalar


#choose path
#options: "circle", "semi_circle", "oval", "eight", "google_earth"
path_name = "circle"

#in case of google_earth specify the .kml
external = 'Map_processing/Maps_kml/extHORTO.kml'
internal = 'Map_processing/Maps_kml/intHORTO.kml'



#Vehicle info
m = 85 #vehicle mass
J=10 #Moment of inertia
mu = 1 #tyre friction coeficient 
pho_air = 1.225 #air density
A0 = 0.5 #frontal area of the car
Cx = 0.5 #Drag coeficient
width = 0.5 #car track width
L = 1 #can wheelbase
Wf=0.4 #position of the center of mass in relation to wheels
h=0.05 #CG height











##################################################################
###                        Path Information                    ###
##################################################################



#create path splines and path spline derivative, assess orientation angles 
#over the sections, define outlines 
spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
    left = choose_path(path_name,external,
internal,N_angle=n_discretization)



#spline points for plotting
spline_points = spline(np.linspace(0,1,num = N_path_points))










##################################################################
###                         Model 1 & 2                        ###
##################################################################


#Define physics over the path. Uncomment the desired Physics model
n_wheels=1 #number of wheels

#Model 1, point
R_t, M_t, C_t, A_t = model1(spline,n_discretization,m,mu)


#Model 2, oriented point with drag
# R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,\
#     pho_air,A0,Cx)



#Comment the models you dont want to compute

#Model abu
# t1_abu=init_optimization_abu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 
"""

#Model bu
t1_bu=init_optimization_bu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 

#Model b
t1_b=init_optimization_b(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)

#Model SOCP abu
t1_SOCP_abu=init_optimization_SOCP_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)

    """
#Model SOCP b
# t1_SOCP_b=init_optimization_SOCP_b(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)




t1_SQP_abu=init_optimization_SQP_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)





"""


##################################################################
###                           Model 3                          ###
##################################################################



#Define physics over the path. 
n_wheels=4 #number of wheels

#Model 3, 4 wheels with drag
R_t, M_t, C_t, A_t = model3(spline,angle,angle_derivative,\
    angle_sec_derivative,n_discretization,m,mu,\
        pho_air,A0,Cx,J,width,L,Wf,n_wheels)



#Comment the models you dont want to compute

#Model abu
t1_abu_3=init_optimization_abu_3(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 


#Model bu
t1_bu_3=init_optimization_bu_3(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 


#Model b
t1_b_3=init_optimization_b_3(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)


#Model SOCP abu
t1_SOCP_abu_3=init_optimization_SOCP_abu_3(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)


#Model SOCP b
t1_SOCP_b_3=init_optimization_SOCP_b_3(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)














##################################################################
###                           Model 4                          ###
##################################################################



#Define physics over the path.
n_wheels=3 #number of wheels


#Model 4, 3 wheels with drag and load transfer

#defines the wheels angles
theta_r,theta_f0,theta_f1 = model4_extra_angles(spline.derivative(),\
    spline.derivative().derivative(),n_discretization,Wf,L,width)

#defines the model's matrices
R_t, M_t, C_t, d_t, A_t = model4(spline,angle,angle_derivative,\
    angle_sec_derivative,theta_r,theta_f0,theta_f1,n_discretization,\
        m,mu,pho_air,A0,Cx,J,width,L,Wf,h,n_wheels)

        

#Comment the models you dont want to compute

#Model abu
t1_abu_4=init_optimization_abu_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 


#Model bu
t1_bu_4=init_optimization_bu_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 


#Model b
t1_b_4=init_optimization_b_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)


#Model SOCP abu
t1_SOCP_abu_4=init_optimization_SOCP_abu_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)


#Model SOCP b
t1_SOCP_b_4=init_optimization_SOCP_b_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)


















##################################################################
###                 Model Performance Comparison              ###
##################################################################

#number of timeit assessments
N_computation_average=1

#Physical model to compute
Physical_model=1

#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]


#Use same order as the models above
#t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]
results = [t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]]


#Call the timeit
#Pay attention to call with the correct physical matrices

d_t = 0 #comment for model 4

computation_time = model_performance(Physical_model,models,results,N_computation_average,
            R_t, M_t,C_t,d_t,A_t,n_discretization,xsi,n_wheels,display=False)

















##################################################################
###                       Comparison Export                   ###
##################################################################

#number of sections to access
discretizations = [10,18,32,56,100]

#number of timeit assessments
N_computation_average=50

#chose the filename
filename = "Comparison/Results/comparison_timeit_model1_eight.csv"

#Physical model to compute
Physical_model=1

#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models_export = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]

export_comparison_to_csv(Physical_model,models_export, discretizations,filename,
    N_computation_average,xsi,spline,m,mu,pho_air,A0,Cx,J,width,L,Wf,h)



#data Dictionary
data = read_csv_to_dict(filename)

#Complexity calculation

complexity = fit_log(data)

title = "Point mass"
model_complexity(models,complexity,title)












##################################################################
###                     Real Path Calculation                  ###
##################################################################


#Only "Time abu" and "Time bu" available
# controlled_path = controlled_path("Time bu",R_t, M_t, C_t, A_t,
#         n_discretization,xsi,n_wheels,spline_points,derivative,N_path_points)


















##################################################################
###                            Plots                           ###
##################################################################


#Uncomment the plots you want

# #solution general model
# t0_abu,t1_abu,forcex0_abu,forcey0_abu,forcex1_abu,forcey1_abu,x0_abu,\
#      decision_variables_abu = init_optimization_abu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)


#  #Test if the circular path velocity is equal to the theoretical
# circular_path_test(derivative,decision_variables_abu[0:n_discretization],n_discretization,m,mu,\
#     pho_air,A0,Cx)





#usar sempre o abu SOCP
t1_SOCP_abu,decision_variables_SOCP_abu = init_optimization_abu_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)

n_wheels = 3



#Animates initial guess vs optimized solution
animation_complete(spline,right,left,spline_points,decision_variables_SOCP_abu.x,\
               t1_SOCP_abu,n_discretization,m,mu,n_wheels)






#Solution comparison plot
# comparison_plot(derivative,R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels)

"""