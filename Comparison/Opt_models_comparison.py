import csv 
import os
from Visualization.print import print_separator
from Simulation.optimization_main import *
from Physics.model1 import model1
import numpy as np
from scipy.stats import linregress




#Exports model results to a CSV file with specified discretizations.
#Input models name (list), discretizations (list), filename (str).

def export_comparison_to_csv(models, discretizations,filename,
                             N_computation_average,xsi,n_wheels,spline,m,mu):
    
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
        R_t, M_t, C_t, A_t = model1(spline,discretizations[i],m,mu)   
        
        #Dictionary with the times traverse
        time_traverse_dict = {}
        
        for name, func in Models_dict.items():
            time_traverse_dict[name] = func(R_t, M_t, C_t, 
            A_t,discretizations[i],xsi,n_wheels,display=False,plot=False)
            
        #save the results 1 list
        Results1.append([time_traverse_dict[name][-1] for name in models])

        #Call the timeit and saves on the second list
        Results2.append(model_performance(models,Results1[i],
            N_computation_average,R_t, M_t, C_t,A_t,discretizations[i],\
                xsi,n_wheels,display=False))
        

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
            for i in range(1, len(row), 3):
                # Get a group of three columns
                triplet = [float(row[i]), float(row[i+1]), float(row[i+2])]  
                values.append(triplet)
                
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
        for i in range(len(results[name])):
            log.append(np.log10(results[name][i][2]))
            n_discretization.append(np.log10(results[name][i][0]))
        
        #fits the straight line    
        slope, intercept, r_value, p_value, std_err = \
            linregress(n_discretization, log)
        
        #store at the dictionary
        dict[name]=[n_discretization,log,[slope, intercept]]
        
            
    return dict