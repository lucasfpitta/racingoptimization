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


n_discretization=20 #number of path sections
N_path_points=1000 #plotting discretization
xsi = 1 #optimization scalar


#choose path
#options: "circle", "semi_circle", "oval", "eight", "google_earth"
path_name = "circle"

#in case of google_earth specify the .kml
external = 'Map_processing/Maps_kml/extHORTO.kml'
internal = 'Map_processing/Maps_kml/intHORTO.kml'



#create path splines and path spline derivative, assess orientation angles 
#over the sections, define outlines 
spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
    left = choose_path(path_name,external,
internal,N_angle=n_discretization)



#spline points for plotting
spline_points = spline(np.linspace(0,1,num = N_path_points))




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
h=0.35 #CG height
 #number of wheels, 1 for model1 and model2, 4 for model3, 3 for model4
n_wheels=3




#Define physics over the path. Uncomment the desired Physics model

#Model 1, point
# R_t, M_t, C_t, A_t = model1(spline,n_discretization,m,mu)

#Model 2, oriented point with drag
# R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,\
#     pho_air,A0,Cx)


#Model 3, 4 wheels with drag
# R_t, M_t, C_t, A_t = model3(spline,angle,angle_derivative,\
#     angle_sec_derivative,n_discretization,m,mu,\
#         pho_air,A0,Cx,J,width,L,Wf,n_wheels)



#Model 4, 3 wheels with drag and load transfer

theta_r,theta_f0,theta_f1 = model4_extra_angles(spline.derivative(),\
    spline.derivative().derivative(),n_discretization,Wf,L,width)


R_t, M_t, C_t, d_t, A_t = model4(spline,angle,angle_derivative,\
    angle_sec_derivative,theta_r,theta_f0,theta_f1,n_discretization,\
        m,mu,pho_air,A0,Cx,J,width,L,Wf,h,n_wheels)






"""

##################################################################
###                           Choose Model                     ###
##################################################################


#Comment the models you dont want to compute

#Model abu
t1_abu=init_optimization_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 


#Model bu
t1_bu=init_optimization_bu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False) 

#Model b
t1_b=init_optimization_b(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)

#Model SOCP abu
t1_SOCP_abu=init_optimization_SOCP_abu(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)

#Model SOCP b
t1_SOCP_b=init_optimization_SOCP_b(
    R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=False)

"""


#Model abu 3
t1_SOCP_abu_4,decision_variables_SOCP_abu_4=init_optimization_SOCP_b_4(
    R_t, M_t, C_t, d_t, A_t,n_discretization,xsi,n_wheels,display=True,plot=True) 

print(t1_SOCP_abu_4)

for i in range(n_wheels):   
    print("lon",decision_variables_SOCP_abu_4[n_discretization+(3*i)*(n_discretization-1):2*n_discretization-1+(3*i)*(n_discretization-1)])
    
for i in range(n_wheels):   
    print("lat",decision_variables_SOCP_abu_4[n_discretization+(3*i+1)*(n_discretization-1):2*n_discretization-1+(3*i+1)*(n_discretization-1)])

for i in range(n_wheels):   
    print("vert", decision_variables_SOCP_abu_4[n_discretization+(3*i+2)*(n_discretization-1):2*n_discretization-1+(3*i+2)*(n_discretization-1)])
    


"""
##################################################################
###                 Model Performance Comparison              ###
##################################################################

#number of timeit assessments
N_computation_average=2


#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]


#Use same order as the models above
#t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]
results = [t1_abu[-1], t1_bu[-1],t1_b[-1],t1_SOCP_abu[-1],t1_SOCP_b[-1]]


#Call the timeit
computation_time = model_performance(models,results,N_computation_average,
            R_t, M_t, C_t,A_t,n_discretization,xsi,n_wheels,display=False)



##################################################################
###                       Comparison Export                   ###
##################################################################

#number of sections to access
discretizations = [5,8]

#number of timeit assessments
N_computation_average=2

#chose the filename
filename = "Comparison/Results/comparison_timeit.csv"

#List to chose the models you do not want to time
#"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"

models_export = ["Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b"]

export_comparison_to_csv(models_export, discretizations,filename,
                         N_computation_average,xsi,n_wheels,spline,m,mu)



#data Dictionary
data = read_csv_to_dict(filename)

#Complexity calculation

complexity = fit_log(data)

model_complexity(models,complexity)


##################################################################
###                     Real Path Calculation                  ###
##################################################################


#Only "Time abu" and "Time bu" available
# controlled_path = controlled_path("Time bu",R_t, M_t, C_t, A_t,
#         n_discretization,xsi,n_wheels,spline_points,derivative,N_path_points)









##################################################################
###                            Plots                           ###
##################################################################
"""

#Uncomment the plots you want

# #solution general model
# t0_abu,t1_abu,forcex0_abu,forcey0_abu,forcex1_abu,forcey1_abu,x0_abu,\
#      decision_variables_abu = init_optimization_abu(
#     R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)


# #Test if the circular path velocity is equal to the theoretical
circular_path_test(derivative,decision_variables_SOCP_abu_4[0:n_discretization],n_discretization,m,mu,\
    pho_air,A0,Cx)

# #Animates initial guess vs optimized solution
# animation_(spline,right,left,spline_points,forcex0_abu,forcey0_abu,\
#             forcex1_abu,forcey1_abu
#                ,t0_abu,t1_abu,n_discretization,m)

#Solution comparison plot
# comparison_plot(derivative,R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels)
