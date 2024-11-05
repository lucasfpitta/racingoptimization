import csv 
import os
from Visualization.print import print_separator
from Simulation.optimization_main import *
from Physics.model1 import model1


#Exports model results to a CSV file with specified discretizations.
#Input models name (list), discretizations (list), filename (str).

def export_comparison_to_csv(models, discretizations,filename,N_computation_average,xsi,spline):
    
    print_separator("Model Comparison")
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    #Defines the time to traverse for each algorithm and discretization
    
     
    
     # Dictionary of Models
    Models_dict = {
        "Time abu": init_optimization_abu,
        "Time bu": init_optimization_bu,
        "Time b": init_optimization_b,
        "Time SOCP abu": init_optimization_SOCP_abu,
        "Time SOCP b": init_optimization_SOCP_b
        }
    
    
    Results1 = []
    Results2 = []
    
    for i in range(len(discretizations)):
        
        #Define physics over the path
        R_t, M_t, C_t, A_t = model1(spline,discretizations[i])   
        
        #Dictionary with the times traverse
        time_traverse_dict = {}
        
        for name, func in Models_dict.items():
            time_traverse_dict[name] = func(R_t, M_t, C_t, 
            A_t,discretizations[i],xsi,display=False,plot=False)
            
        #save the results 1 list
        Results1.append([time_traverse_dict[name][-1] for name in models])

        #Call the timeit and saves on the second list
        Results2.append(model_performance(models,Results1[i],
            N_computation_average,R_t, M_t, C_t,A_t,discretizations[i],xsi,
            display=False))
        

    # Create headers with dynamic discretization triplets
    headers = ["Model Name"]
    
    for d in discretizations:
        headers.extend([f"N = {d}", f"Time Traverse N = {d}", 
                        f"Time Compute N= {d}"])

    # Write to CSV
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for i in range(len(models)):
            row = [models[i]]
            for d in range(len(discretizations)):
                row.extend([discretizations[d], Results1[d][i], Results2[d][i]])
            writer.writerow(row)
            
    print()
    print(f"Results exported to {filename}.")
    
    return