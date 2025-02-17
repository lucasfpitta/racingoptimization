import tkinter as tk
import numpy as np
import sys
import ast

from tkinter import ttk

from Map_processing.choose_path import choose_path
from Physics.model1 import model1
from Physics.model2 import model2
from Physics.model3 import model3
from Physics.model4 import model4
from splines.splines import model4_extra_angles
from Simulation.optimization_main import *
from Comparison.Opt_models_comparison import *
from Visualization.plots import *




##################################################################
###                         Test Functions                     ###
##################################################################

config = None


def test_model_1():
    
    """"
    Path information
    """
    
    #create path splines and path spline derivative, assess orientation angles 
    #over the sections, define outlines 
    spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
        left = choose_path(config.path_name,config.external,config.internal,
        N_angle=config.n_discretization)


    #spline points for plotting
    spline_points = spline(np.linspace(0,1,num = config.N_path_points))
    
    
    
    """"
    Optimization
    """
    
    #Define physics over the path. Uncomment the desired Physics model
    n_wheels=1 #number of wheels

    #Model 1, point
    R_t, M_t, C_t, A_t = model1(spline,config.n_discretization,config.m,config.mu)

    #Comment the models you dont want to compute

    #Model abu
    t1_abu=init_optimization_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 


    #Model bu
    t1_bu=init_optimization_bu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 

    #Model b
    t1_b=init_optimization_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP abu
    t1_SOCP_abu=init_optimization_SOCP_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP b
    t1_SOCP_b=init_optimization_SOCP_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SQP trust-constr abu
    t1_SQP_abu=init_optimization_SQP_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)

    #Model SQP trust-constr b
    t1_SQP_b=init_optimization_SQP_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


def test_model_2():
    
    
    """"
    Path information
    """
    
    #create path splines and path spline derivative, assess orientation angles 
    #over the sections, define outlines 
    spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
        left = choose_path(config.path_name,config.external,config.internal,
        N_angle=config.n_discretization)


    #spline points for plotting
    spline_points = spline(np.linspace(0,1,num = config.N_path_points))
    
    
    
    """"
    Optimization
    """
    
    #Define physics over the path. Uncomment the desired Physics model
    n_wheels=1 #number of wheels

    #Model 2, oriented point with drag
    R_t, M_t, C_t, A_t = model2(spline,angle,config.n_discretization,config.m,config.mu,\
    config.pho_air,config.A0,config.Cx)

    #Comment the models you dont want to compute

    #Model abu
    t2_abu=init_optimization_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 


    #Model bu
    t2_bu=init_optimization_bu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 

    #Model b
    t2_b=init_optimization_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP abu
    t2_SOCP_abu=init_optimization_SOCP_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP b
    t2_SOCP_b=init_optimization_SOCP_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SQP trust-constr abu
    t2_SQP_abu=init_optimization_SQP_abu(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)

    #Model SQP trust-constr b
    t2_SQP_b=init_optimization_SQP_b(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)







def test_model_3():
    
    
    """"
    Path information
    """
    
    #create path splines and path spline derivative, assess orientation angles 
    #over the sections, define outlines 
    spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
        left = choose_path(config.path_name,config.external,config.internal,
        N_angle=config.n_discretization)


    #spline points for plotting
    spline_points = spline(np.linspace(0,1,num = config.N_path_points))
    
    
    
    """"
    Optimization
    """
    
    
    
    
    #Define physics over the path. 
    n_wheels=4 #number of wheels

    #Model 3, 4 wheels with drag
    R_t, M_t, C_t, A_t = model3(spline,angle,angle_derivative,\
        angle_sec_derivative,config.n_discretization,config.m,config.mu,\
        config.pho_air,config.A0,config.Cx,config.J,config.width,config.L,config.Wf,n_wheels)



    #Comment the models you dont want to compute

    #Model abu
    t1_abu_3=init_optimization_abu_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 


    #Model bu
    t1_bu_3=init_optimization_bu_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 


    #Model b
    t1_b_3=init_optimization_b_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP abu
    t1_SOCP_abu_3=init_optimization_SOCP_abu_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SOCP b
    t1_SOCP_b_3=init_optimization_SOCP_b_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)



    #Model SQP abu
    t1_SQP_abu_3=init_optimization_SQP_abu_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SQP b
    t1_SQP_b_3=init_optimization_SQP_b_3(
        R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)



def test_model_4():
    
    
    """"
    Path information
    """
    
    #create path splines and path spline derivative, assess orientation angles 
    #over the sections, define outlines 
    spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
        left = choose_path(config.path_name,config.external,config.internal,
        N_angle=config.n_discretization)


    #spline points for plotting
    spline_points = spline(np.linspace(0,1,num = config.N_path_points))
    
    
    
    """"
    Optimization
    """
    
    
    
    
    #Define physics over the path.
    n_wheels=3 #number of wheels


    #Model 4, 3 wheels with drag and load transfer

    #defines the wheels angles
    theta_r,theta_f0,theta_f1 = model4_extra_angles(spline.derivative(),\
        spline.derivative().derivative(),config.n_discretization,config.Wf,config.L,config.width)

    #defines the model's matrices
    R_t, M_t, C_t, d_t, A_t = model4(spline,angle,angle_derivative,\
        angle_sec_derivative,theta_r,theta_f0,theta_f1,config.n_discretization,\
    config.m,config.mu,config.pho_air,config.A0,config.Cx,config.J,config.width,\
        config.L,config.Wf,config.h,n_wheels)

            

    #Comment the models you dont want to compute

    #Model abu
    t1_abu_4=init_optimization_abu_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 



    #Model bu
    t1_bu_4=init_optimization_bu_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False) 



    #Model b
    t1_b_4=init_optimization_b_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)



    #Model SOCP abu
    t1_SOCP_abu_4=init_optimization_SOCP_abu_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)



    #Model SOCP b
    t1_SOCP_b_4=init_optimization_SOCP_b_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)


    #Model SQP abu
    t1_SQP_abu_4=init_optimization_SQP_abu_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)



    #Model SQP b
    t1_SQP_b_4=init_optimization_SQP_b_4(
        R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=False)




def model_comparison_export():
    
    """"
    Path information
    """
    
    #create path splines and path spline derivative, assess orientation angles 
    #over the sections, define outlines 
    spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
        left = choose_path(config.path_name,config.external,config.internal,
        N_angle=config.n_discretization)
    
    export_comparison_to_csv(config.Physical_model,config.models_export,config.discretizations,\
    config.filename,config.N_computation_average,config.xsi,spline,config.m,config.mu,\
    config.pho_air,config.A0,config.Cx,config.J,config.width,config.L,config.Wf,config.h)
    
    #data Dictionary
    data = read_csv_to_dict(config.filename)

    #Complexity calculation

    complexity = fit_log(data)


    #"Time abu","Time bu","Time b","Time SOCP abu","Time SOCP b","Time SQP abu","Time SQP b"
    title = "Models Complexity Comparison"
    model_complexity(config.models_export,complexity,title)
    
    
    
    
    
    

def plots():
    
    
    
    if config.Physical_model==1 or config.Physical_model ==2:
        #create path splines and path spline derivative, assess orientation angles 
        #over the sections, define outlines 
        spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
            left = choose_path(config.path_name,config.external,
        config.internal,N_angle=config.n_discretization)
            
         #Define physics over the path. Uncomment the desired Physics model
        n_wheels=1 #number of wheels

        if config.Physical_model==1:
            #Model 1, point
            R_t, M_t, C_t, A_t = model1(spline,config.n_discretization,config.m,config.mu)
        elif config.Physical_model==2:
            #Model 2, oriented point with drag
            R_t, M_t, C_t, A_t = model2(spline,angle,config.n_discretization,config.m,config.mu,\
                config.pho_air,config.A0,config.Cx)


        #spline points for plotting
        spline_points = spline(np.linspace(0,1,num = config.N_path_points))
        
        t1_SOCP_abu,decision_variables_SOCP_abu = init_optimization_SOCP_abu(
            R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=False,plot=True)




    elif config.Physical_model==3:
        #create path splines and path spline derivative, assess orientation angles 
        #over the sections, define outlines 
        spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
            left = choose_path(config.path_name,config.external,config.internal,
            N_angle=config.n_discretization)

        #spline points for plotting
        spline_points = spline(np.linspace(0,1,num = config.N_path_points))
        
        #Define physics over the path. 
        n_wheels=4 #number of wheels

        #Model 3, 4 wheels with drag
        R_t, M_t, C_t, A_t = model3(spline,angle,angle_derivative,\
            angle_sec_derivative,config.n_discretization,config.m,config.mu,\
            config.pho_air,config.A0,config.Cx,config.J,config.width,config.L,config.Wf,n_wheels)

        #Model SOCP abu
        t1_SOCP_abu,decision_variables_SOCP_abu=init_optimization_SOCP_abu_3(
            R_t, M_t, C_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=True)

    
    
    elif config.Physical_model==4:
        
        #create path splines and path spline derivative, assess orientation angles 
        #over the sections, define outlines 
        spline, derivative, angle, angle_derivative,angle_sec_derivative, right,\
            left = choose_path(config.path_name,config.external,config.internal,
            N_angle=config.n_discretization)


        #spline points for plotting
        spline_points = spline(np.linspace(0,1,num = config.N_path_points))
 
        #Define physics over the path.
        n_wheels=3 #number of wheels

        #defines the wheels angles
        theta_r,theta_f0,theta_f1 = model4_extra_angles(spline.derivative(),\
            spline.derivative().derivative(),config.n_discretization,config.Wf,config.L,config.width)

        #defines the model's matrices
        R_t, M_t, C_t, d_t, A_t = model4(spline,angle,angle_derivative,\
            angle_sec_derivative,theta_r,theta_f0,theta_f1,config.n_discretization,\
        config.m,config.mu,config.pho_air,config.A0,config.Cx,config.J,config.width,\
            config.L,config.Wf,config.h,n_wheels)


        #Model SOCP abu
        t1_SOCP_abu,decision_variables_SOCP_abu=init_optimization_SOCP_abu_4(
            R_t, M_t, C_t, d_t, A_t,config.n_discretization,config.xsi,n_wheels,display=True,plot=True)

    else:
        print("Wrong Physical Model Number")
        
    #Animates initial guess vs optimized solution
    animation_complete(spline,right,left,spline_points,decision_variables_SOCP_abu,\
                    t1_SOCP_abu,config.n_discretization,config.m,config.mu,n_wheels)
    
    #compares local max velocity and optimize velocity
    local_max_v(derivative,decision_variables_SOCP_abu[config.n_discretization-1:2*config.n_discretization-1]\
    ,config.n_discretization,config.m,config.mu,config.pho_air,config.A0,config.Cx)








##################################################################
###                       Graphical Interface                  ###
##################################################################


def screen():
    
     #Custom ScrollableFrame Class Code
    class ScrollableFrame(ttk.Frame):
        def __init__(self, container, width=500, height=150, *args, **kwargs):
            super().__init__(container, *args, **kwargs)
            self.canvas = tk.Canvas(self, width=width, height=height)
            self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
            self.scrollable_frame = ttk.Frame(self.canvas)

            # Update the scrollregion when the size of the scrollable frame changes
            self.scrollable_frame.bind(
                "<Configure>",
                lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
            )
            # Add the scrollable frame to the canvas
            self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
            self.canvas.configure(yscrollcommand=self.scrollbar.set)

            # Pack canvas and scrollbar
            self.canvas.pack(side="left", fill="both", expand=True)
            self.scrollbar.pack(side="right", fill="y")

 
    class RedirectText:
        def __init__(self, text_widget):
            self.text_widget = text_widget

        def write(self, string):
            self.text_widget.insert(tk.END, string)
            self.text_widget.see(tk.END)

        def flush(self):
            pass
 
    #SUbparameters for google earth map
    subparam_entries_maps = {}
    def google_earth_map(*args):
   
        #Clears and repopulates the subparameters frame based on the selected function.

        global subparam_entries_maps
        # Clear previous entries
        subparam_entries_maps = {}
        for widget in scroll_file.scrollable_frame.winfo_children():
            widget.destroy()

        map_selection = path_var.get()
        
        if map_selection == "google_earth":

            tk.Label(scroll_file.scrollable_frame, text="external (kml file):").grid(row=1, column=0, sticky="w")
            entry_external = ttk.Combobox(scroll_file.scrollable_frame, values=["Map_processing/Maps_kml/extHORTO.kml"], width=40)
            entry_external.set("Map_processing/Maps_kml/extHORTO.kml")
            entry_external.grid(row=1, column=1, padx=5, pady=2)
            subparam_entries_maps['external'] = entry_external

            tk.Label(scroll_file.scrollable_frame, text="internal (kml file):").grid(row=2, column=0, sticky="w")
            entry_internal = ttk.Combobox(scroll_file.scrollable_frame, values=["Map_processing/Maps_kml/intHORTO.kml"], width=40)
            entry_internal.set("Map_processing/Maps_kml/intHORTO.kml")
            entry_internal.grid(row=2, column=1, padx=5, pady=2)
            subparam_entries_maps['internal'] = entry_internal
            
        return subparam_entries_maps


    subparam_entries_func = {}       
    #Gets subparmeters for functions
    def func_info(*args):
   
        global subparam_entries_func
        # Clear previous entries
        subparam_entries_func = {}
        for widget in scroll_func_info.scrollable_frame.winfo_children():
            widget.destroy()

        func_selection = func_var.get()
        
        model_list = ["model_1","model_2","model_3","model_4"]
        if func_selection in model_list:

            tk.Label(scroll_func_info.scrollable_frame, text="# Path sections (int):").grid(row=0, column=0, sticky="w")
            entry_n_discretization = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_n_discretization.insert(0, "33")
            entry_n_discretization.grid(row=1, column=1, padx=5, pady=2)
            subparam_entries_func['n_discretization'] = entry_n_discretization 
            
        if func_selection == "plots":

            tk.Label(scroll_func_info.scrollable_frame, text="# Path sections (int):").grid(row=0, column=0, sticky="w")
            entry_n_discretization = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_n_discretization.insert(0, "200")
            entry_n_discretization.grid(row=0, column=1, padx=5, pady=2)
            subparam_entries_func['n_discretization'] = entry_n_discretization 
            
            tk.Label(scroll_func_info.scrollable_frame, text="# Physical model (int):").grid(row=1, column=0, sticky="w")
            entry_Physical_model = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_Physical_model.insert(0, "3")
            entry_Physical_model.grid(row=1, column=1, padx=5, pady=2)
            subparam_entries_func['Physical_model'] = entry_Physical_model 
            
            
        elif func_selection == "comparison_export":

            tk.Label(scroll_func_info.scrollable_frame, text="# Path sections (list):").grid(row=0, column=0, sticky="w")
            entry_n_discretization = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_n_discretization.insert(0, "[2,3]")
            entry_n_discretization.grid(row=0, column=1, padx=5, pady=2)
            subparam_entries_func['discretizations'] = entry_n_discretization 
            
            tk.Label(scroll_func_info.scrollable_frame, text="# Comput assess (int):").grid(row=1, column=0, sticky="w")
            entry_n_avg = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_n_avg.insert(0, "10")
            entry_n_avg.grid(row=1, column=1, padx=5, pady=2)
            subparam_entries_func['n_avg'] = entry_n_avg
            
            tk.Label(scroll_func_info.scrollable_frame, text="# Physical model (int):").grid(row=2, column=0, sticky="w")
            entry_Physical_model = tk.Entry(scroll_func_info.scrollable_frame, width=10)
            entry_Physical_model.insert(0, "3")
            entry_Physical_model.grid(row=2, column=1, padx=5, pady=2)
            subparam_entries_func['Physical_model'] = entry_Physical_model 
            
        return subparam_entries_func  
        
        




    # --- Function to retrieve parameters ---
    def run_parameters():
        global subparam_entries_maps
        global subparam_entries_func
        
        try:
            # Optimization variables:
            config.N_path_points = int(entry_N_path_points.get())
            config.xsi = float(entry_xsi.get())
            
            # Path selection:
            config.path_name = path_var.get()
            if config.path_name == "google_earth":
                subparam_entries_maps = google_earth_map()
                config.external = subparam_entries_maps['external'].get()
                config.internal = subparam_entries_maps['internal'].get()   

            # Vehicle info:
            config.m = float(entry_m.get())
            config.J = float(entry_J.get())
            config.mu = float(entry_mu.get())
            config.pho_air = float(entry_pho_air.get())
            config.A0 = float(entry_A0.get())
            config.Cx = float(entry_Cx.get())
            config.width = float(entry_width.get())
            config.L = float(entry_L.get())
            config.Wf = float(entry_Wf.get())
            config.h = float(entry_h.get())
            
            function = func_var.get()
            
            if function == "model_1":
                config.n_discretization = int(subparam_entries_func['n_discretization'].get())
                test_model_1()
                
            if function == "model_2":
                config.n_discretization = int(subparam_entries_func['n_discretization'].get())
                test_model_2()
                
            if function == "model_3":
                config.n_discretization = int(subparam_entries_func['n_discretization'].get())
                test_model_3()
                
            if function == "model_4":
                config.n_discretization = int(subparam_entries_func['n_discretization'].get())
                test_model_4()
                
            if function == "plots":
                config.n_discretization = int(subparam_entries_func['n_discretization'].get())
                config.Physical_model = int(subparam_entries_func['Physical_model'].get())
                plots()
                
            if function == "comparison_export":
                config.discretizations = ast.literal_eval(subparam_entries_func['discretizations'].get())
                config.N_computation_average = int(subparam_entries_func['n_avg'].get())
                config.Physical_model = int(subparam_entries_func['Physical_model'].get())
                model_comparison_export()
        
        
            #result_label.config(text="Check your terminal")
        except ValueError:
            print("Error: Please check that all inputs are valid numbers.")

    # --- Main Window Setup ---
    root = tk.Tk()
    root.title("Parameter Input")
    
    
    
    
    
    
    
    # --- Vehicle Info Frame  ---
    frame_vehicle = ttk.LabelFrame(root, text="Vehicle Info", padding=10)
    frame_vehicle.grid(row=0, column=0, rowspan=3, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Vehicle Info frame
    scroll_vehicle = ScrollableFrame(frame_vehicle, width=300, height=250)
    scroll_vehicle.pack(fill="both", expand=True)

    tk.Label(scroll_vehicle.scrollable_frame, text="mass (kg):").grid(row=0, column=0, sticky="w")
    entry_m = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_m.insert(0, "85")
    entry_m.grid(row=0, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="Moment of inertia (kg mˆ2):").grid(row=1, column=0, sticky="w")
    entry_J = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_J.insert(0, "10")
    entry_J.grid(row=1, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="tyre friction coef:").grid(row=2, column=0, sticky="w")
    entry_mu = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_mu.insert(0, "1")
    entry_mu.grid(row=2, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="air density (kg/mˆ3):").grid(row=3, column=0, sticky="w")
    entry_pho_air = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_pho_air.insert(0, "1.225")
    entry_pho_air.grid(row=3, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="Frontal area (mˆ2):").grid(row=4, column=0, sticky="w")
    entry_A0 = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_A0.insert(0, "0.5")
    entry_A0.grid(row=4, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="Drag coef.:").grid(row=5, column=0, sticky="w")
    entry_Cx = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_Cx.insert(0, "0.5")
    entry_Cx.grid(row=5, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="Track width (m):").grid(row=6, column=0, sticky="w")
    entry_width = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_width.insert(0, "0.5")
    entry_width.grid(row=6, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="Wheelbase (m):").grid(row=7, column=0, sticky="w")
    entry_L = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_L.insert(0, "1")
    entry_L.grid(row=7, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="CG position ([0,1]):").grid(row=8, column=0, sticky="w")
    entry_Wf = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_Wf.insert(0, "0.4")
    entry_Wf.grid(row=8, column=1, padx=5, pady=2)

    tk.Label(scroll_vehicle.scrollable_frame, text="CG height (m):").grid(row=9, column=0, sticky="w")
    entry_h = tk.Entry(scroll_vehicle.scrollable_frame, width=10)
    entry_h.insert(0, "0.35")
    entry_h.grid(row=9, column=1, padx=5, pady=2)
    
    
    

    # --- Optimization Variables Frame (Left Bottom) ---
    frame_opt = ttk.LabelFrame(root, text="Optimization Variables", padding=10)
    frame_opt.grid(row=3, column=0, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Optimization Variables frame
    scroll_opt = ScrollableFrame(frame_opt, width=300, height=100)
    scroll_opt.pack(fill="both", expand=True)

    tk.Label(scroll_opt.scrollable_frame, text="# Plot discretizations (int):").grid(row=0, column=0, sticky="w")
    entry_N_path_points = tk.Entry(scroll_opt.scrollable_frame, width=10)
    entry_N_path_points.insert(0, "1000")
    entry_N_path_points.grid(row=0, column=1, padx=5, pady=2)

    tk.Label(scroll_opt.scrollable_frame, text="xsi (float [0,1]):").grid(row=1, column=0, sticky="w")
    entry_xsi = tk.Entry(scroll_opt.scrollable_frame, width=10)
    entry_xsi.insert(0, "1")
    entry_xsi.grid(row=1, column=1, padx=5, pady=2)









    # --- Path Selection Frame ---
    frame_path = ttk.LabelFrame(root, text="Path Selection", padding=10)
    frame_path.grid(row=0, column=1, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Path Selection frame
    scroll_path = ScrollableFrame(frame_path, width=300, height=50)
    scroll_path.pack(fill="both", expand=True)
    
    
    path_var = tk.StringVar(value="eight")              

    tk.Label(scroll_path.scrollable_frame, text="Path name:").grid(row=0, column=0, sticky="w")
    combo_path_name = ttk.Combobox(scroll_path.scrollable_frame, textvariable=path_var,values=["circle", "semi_circle", "oval", "eight", "google_earth"], width=15)
    combo_path_name.set("eight")
    combo_path_name.grid(row=0, column=1, padx=5, pady=2)
    
    
    
    
    
    
    # --- Path Selection File Frame   ---
    frame_file = ttk.LabelFrame(root, text="File Selection", padding=10)
    frame_file.grid(row=1, column=1, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Path Selection frame
    scroll_file = ScrollableFrame(frame_file, width=300, height=50)
    scroll_file.pack(fill="both", expand=True)              

    # Update subparameters whenever the function selection changes
    path_var.trace("w", google_earth_map)






# --- Function Selection Frame  ---
    frame_func = ttk.LabelFrame(root, text="Function Selection", padding=10)
    frame_func.grid(row=2, column=1, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Path Selection frame
    scroll_func = ScrollableFrame(frame_func, width=300, height=50)
    scroll_func.pack(fill="both", expand=True)
    
    
    func_var = tk.StringVar(value="model_1")              

    tk.Label(scroll_func.scrollable_frame, text="Function name:").grid(row=0, column=0, sticky="w")
    combo_func_name = ttk.Combobox(scroll_func.scrollable_frame, textvariable=func_var,values=["model_1","model_2","model_3","model_4","comparison_export","plots"], width=15)
    combo_func_name.set("model_1")
    combo_func_name.grid(row=0, column=1, padx=5, pady=2)
    




    # --- Function additional info Frame   ---
    frame_func_info = ttk.LabelFrame(root, text="Function Variables", padding=10)
    frame_func_info.grid(row=3, column=1, padx=10, pady=5, sticky="nsew")

    # Create a scrollable area within the Path Selection frame
    scroll_func_info = ScrollableFrame(frame_func_info, width=300, height=50)
    scroll_func_info.pack(fill="both", expand=True)              

    # Update subparameters whenever the function selection changes
    func_var.trace("w", func_info)





    # Create the scrollable frame inside root
    scroll_result_label = ScrollableFrame(root, width=500, height=1000)
    scroll_result_label.grid(row=0, column=2, columnspan=4, padx=10, pady=10, sticky="nsew")

    # Add a Text widget inside the scrollable frame
    result_text = tk.Text(scroll_result_label.scrollable_frame, wrap="word", width=100, height=100)
    result_text.pack(fill="both", expand=True)
  
    
    # Redirect stdout to our Text widget
    sys.stdout = RedirectText(result_text)
    
    print("Results")

    # --- Submit Button and Results ---
    submit_button = tk.Button(root, text="Submit", command=run_parameters)
    submit_button.grid(row=4, column=1, padx=10, pady=10, sticky="w")

    
    

    # Allow resizing
    root.columnconfigure(0, weight=1)
    root.columnconfigure(1, weight=1)
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)

    root.mainloop()
