import csv 
import os
from Visualization.print import print_separator
from Simulation.optimization_main import *
from Physics.model1 import model1
from Physics.model2 import model2
from Physics.model3 import model3
from Physics.model4 import model4
import numpy as np
from scipy.stats import linregress
from splines.splines import model4_extra_angles
from splines.splines import find_angle




#Exports model results to a CSV file with specified discretizations.
#Input models name (list), discretizations (list), filename (str).

def export_comparison_to_csv(Physical_model,models, discretizations,filename,
    N_computation_average,xsi,spline,m,mu,pho_air,A0,Cx,J,width,L,Wf,h):
    
    print_separator("Model Comparison")
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    #Defines the time to traverse for each algorithm and discretization
    
     
    
    # Dictionary of Models
     
    if Physical_model == 1 or Physical_model == 2:
        Models_dict = {
               "Time abu": init_optimization_abu,
               "Time bu": init_optimization_bu,
               "Time b": init_optimization_b,
               "Time SOCP abu": init_optimization_SOCP_abu,
               "Time SOCP b": init_optimization_SOCP_b,
               "Time SQP abu": init_optimization_SQP_abu,
               "Time SQP b": init_optimization_SQP_b
               }
        
    elif Physical_model == 3:
        Models_dict = {
               "Time abu": init_optimization_abu_3,
               "Time bu": init_optimization_bu_3,
               "Time b": init_optimization_b_3,
               "Time SOCP abu": init_optimization_SOCP_abu_3,
               "Time SOCP b": init_optimization_SOCP_b_3,
               "Time SQP abu": init_optimization_SQP_abu_3,
               "Time SQP b": init_optimization_SQP_b_3
               }
        
    elif Physical_model == 4:
        Models_dict = {
               #"Time abu": init_optimization_abu_4,
               #"Time bu": init_optimization_bu_4,
               #"Time b": init_optimization_b_4,
               "Time SOCP abu": init_optimization_SOCP_abu_4,
               "Time SOCP b": init_optimization_SOCP_b_4,
               #"Time SQP abu": init_optimization_SQP_abu_4,
               #"Time SQP b": init_optimization_SQP_b_4
               }
        
    else:
        print("Wrong Physical Model")
    
    
    Results1 = []
    Results2 = []
    Results3 = []
    
    for i in range(len(discretizations)):

        angle, angle_derivative, angle_sec_derivative = \
            find_angle(spline,discretizations[i])
        
        #Dictionary with the times traverse
        time_traverse_dict = {}
        
        
        if Physical_model == 1:
            #Define physics over the path
            R_t, M_t, C_t, A_t = model1(spline,discretizations[i],m,mu)   
            
            for name, func in Models_dict.items():
                time_traverse_dict[name] = func(R_t, M_t, C_t, 
                A_t,discretizations[i],xsi,n_wheels=1,display=False,plot=False)

            Results1.append([time_traverse_dict[name][-1] for name in models])
            #Call the timeit and saves on the second list
            d_t=0
            mean,std = model_performance(Physical_model,models,Results1[i],\
                N_computation_average,R_t, M_t, C_t,d_t,A_t,discretizations[i],\
                xsi,n_wheels=1,display=False)
                
                
                
        if Physical_model == 2:
            #Define physics over the path
            R_t, M_t, C_t, A_t = model2(spline,angle,discretizations[i],m,mu,\
                pho_air,A0,Cx)   
            
            for name, func in Models_dict.items():
                time_traverse_dict[name] = func(R_t, M_t, C_t, 
                A_t,discretizations[i],xsi,n_wheels=1,display=False,plot=False)
            
            Results1.append([time_traverse_dict[name][-1] for name in models])
            #Call the timeit and saves on the second list
            d_t=0
            mean,std = model_performance(Physical_model,models,Results1[i],
            N_computation_average,R_t, M_t, C_t,d_t,A_t,discretizations[i],\
                xsi,n_wheels=1,display=False)
                
                
        if Physical_model == 3:
            #Define physics over the path
            R_t, M_t, C_t, A_t=model3(spline,angle,angle_derivative,\
                angle_sec_derivative,discretizations[i],m,mu,\
                pho_air,A0,Cx,J,width,L,Wf,n_wheels=4)   
            
            for name, func in Models_dict.items():
                time_traverse_dict[name] = func(R_t, M_t, C_t, A_t,\
                discretizations[i],xsi,n_wheels=4,display=False,plot=False)
            
            Results1.append([time_traverse_dict[name][-1] for name in models])
            #Call the timeit and saves on the second list
            d_t=0
            mean,std = model_performance(Physical_model,models,Results1[i],
            N_computation_average,R_t, M_t, C_t,d_t,A_t,discretizations[i],\
                xsi,n_wheels=4,display=False)
                
        if Physical_model == 4:
            #defines the wheels angles
            theta_r,theta_f0,theta_f1 = model4_extra_angles(spline.derivative(),\
                spline.derivative().derivative(),discretizations[i],Wf,L,width)

            #defines the model's matrices
            R_t, M_t, C_t, d_t, A_t = model4(spline,angle,angle_derivative,\
                angle_sec_derivative,theta_r,theta_f0,theta_f1,discretizations[i],\
                    m,mu,pho_air,A0,Cx,J,width,L,Wf,h,n_wheels=3)
            
            for name, func in Models_dict.items():
                time_traverse_dict[name] = func(R_t, M_t, C_t,d_t, 
                A_t,discretizations[i],xsi,n_wheels=3,display=False,plot=False)
            
            Results1.append([time_traverse_dict[name][-1] for name in models])
            #Call the timeit and saves on the second list
            mean,std = model_performance(Physical_model,models,Results1[i],
            N_computation_average,R_t, M_t, C_t,d_t,A_t,discretizations[i],\
                xsi,n_wheels=3,display=False)
            
            
        
            
        #save the results 1 list
        Results2.append(mean)
        Results3.append(std)
        

    # Create headers with dynamic discretization triplets
    headers = ["Model Name"]
    
    for d in discretizations:
        headers.extend([f"N = {d}", f"Time Traverse N = {d}", 
                        f"Time Compute N= {d}", f"Std Compute N= {d}"])

    # Write to CSV
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for i in range(len(models)):
            row = [models[i]]
            for d in range(len(discretizations)):
                row.extend([discretizations[d],Results1[d][i],Results2[d][i],Results3[d][i]])
            writer.writerow(row)
            
    print()
    print(f"Results exported to {filename}.")
    
    return












#Reads a CSV file and stores it in a dictionary.
#Input filename The name of the CSV file to read.
#A dictionary where the first column is the key, and each key has a list 
# of lists containing groups of three columns from the rest of the row.
def read_csv_to_dict(filename):

    data_dict = {}
    
    with open(filename, mode="r") as file:
        csv_reader = csv.reader(file)
        
        # Skip the first line (header)
        next(csv_reader)
        
        # Process each row
        for row in csv_reader:
            # Use the first column as the dictionary key
            key = row[0]
            
            # Group remaining columns in triplets
            values = []
            for i in range(1, len(row), 4):
                # Get a group of three columns
                quadruplet = [float(row[i]), float(row[i+1]), float(row[i+2]), float(row[i+3])]  
                values.append(quadruplet)
                
            # Add to dictionary
            data_dict[key] = values
    
    return data_dict
















#takes the logarithm on the time to compute and fits a straight line
#input results dictionary (from csv reader)
#Output all dict model name, log points, 2 line coef

def fit_log(results):
    
    dict = {}
    
    #for each algorithm
    for name in results:

        #Take the log of time to compute
        log = []

        #Take the n_discretization
        n_discretization = []

        #take the confidense interval on the log
        lb = []
        ub = []
        std = []

        for i in range(len(results[name])):
            log.append(np.log10(results[name][i][2]))
            n_discretization.append(np.log10(results[name][i][0]))
            lb.append(np.log10(results[name][i][2]-1.96*results[name][i][3]))
            ub.append(np.log10(results[name][i][2]+1.96*results[name][i][3]))
            std.append(results[name][i][3])
            
        
        #fits the straight line    
        slope, intercept, r_value, p_value, std_err = \
            linregress(n_discretization, log)
        
        #store at the dictionary
        dict[name]=[n_discretization,log,[slope, intercept],lb,ub,std]
        
            
    return dict