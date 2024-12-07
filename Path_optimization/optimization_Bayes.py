import numpy as np
from Map_processing.mapreader import get_outline
from splines.splines import path_info
from Physics.model2 import model2
from Simulation.optimization_main import init_optimization_SOCP_b
import matplotlib.pyplot as plt
from Visualization.plots import animation_complete
from scipy.optimize import differential_evolution



from skopt import gp_minimize
from skopt.space import Real




class LivePlotter:
    def __init__(self):
        self.iteration_data = []
        self.cost_data = []
        self.best_cost = []
        
        # Set up the plot
        plt.ion()  # Interactive mode
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], label="Cost", color='blue', linewidth=2)
        self.line2, = self.ax.plot([], [], label="Best Cost", color='red', linewidth=2)
        self.ax.set_title("Optimization Convergence")
        self.ax.set_xlabel("Iteration")
        self.ax.set_ylabel("Cost")
        self.ax.legend()
        self.ax.grid()

    def update(self, iteration, cost):
        # Update the data
        self.iteration_data.append(iteration)
        self.cost_data.append(cost)
        
        if len(self.best_cost)==0 or cost<self.best_cost[-1]:
            self.best_cost.append(cost)
        else:
            self.best_cost.append(self.best_cost[-1])

        # Update the plot
        self.line.set_data(self.iteration_data, self.cost_data)
        self.line2.set_data(self.iteration_data, self.best_cost)
        self.ax.relim()  # Recalculate limits
        self.ax.autoscale_view()  # Rescale view
        plt.pause(0.01)  # Pause to refresh the plot

    def finalize(self):
        plt.ioff()
        plt.show()









def init_path_optimization(right,left,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels):

    def path_optimization(alfas, plotter=None):
        alfa_v = np.concatenate((alfas,[alfas[0]]))
        spline, derivative, angle, angle_derivative, angle_sec_derivative = \
        path_info(left, right, alfa_v,N_angle)
        R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
        t1_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=False)
        cost = t1_SOCP_b[-1]

        print("current cost: ", cost)
        if plotter is not None:
            iteration = len(plotter.iteration_data)  # Determine current iteration
            plotter.update(iteration, cost)

        return cost
    return path_optimization


def trajectory_optimization(external,internal,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels):
    #defines the outline
    right, left,  = get_outline(external,internal)

    # Define problem dimensions and bounds
    dimensions = len(right[0])-1  # Number of control_points
    print(f"The number of control points is {dimensions}")
   # Define the search space (8 variables between 0 and 1)
    space = [Real(0, 1, name=f"x{i}") for i in range(dimensions)]

    objective = init_path_optimization(right,left,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels)
    # Create an instance of the tracker
    #tracker = BestResultTracker()

    plotter = LivePlotter()
    def wrapped_objective(x):
        return objective(x, plotter)

    # Run Bayesian optimization
    result_b = gp_minimize(
        func=wrapped_objective,  # Objective function to minimize
        dimensions=space,         # Search space
        acq_func="EI",            # Acquisition function: Expected Improvement
        n_calls=200,               # Number of evaluations
        random_state=42           # Seed for reproducibility
    )

    # Print the results
    print("Best score Bayes:", result_b.fun)
    print("Best solution Bayes:", result_b.x)
    np.savetxt('output_path_bayes.txt', result_b.x, fmt='%.3f', delimiter=',', header='Best PATH Bayes', comments='')

    new_bounds=[]
    for i in range(dimensions):
        new_bounds.append((max(result_b.x[i]-0.05,0),min(result_b.x[i]+0.05,1)))

     # Expand the Bayesian result to form a valid population
    population_size = 10  

    # Run Differential Evolution
    result_de = differential_evolution(
        func=wrapped_objective,
        bounds=new_bounds,
        #init=initial_population,  # Use the valid population format
        strategy="best2bin",
        mutation=(0.5, 0.9),  # Mutation range for more diverse steps
        recombination=0.7,    # Default recombination
        maxiter=100,
        popsize=population_size,  # Ensure consistency with init size
        tol=1e-6,
        seed=42
    )

    # Print the results
    print("Best score Differential Evolution:", result_de.fun)
    print("Best solution Differential Evolution:", result_de.x)

    plotter.finalize()
    alfa_v = np.concatenate((result_de.x,[result_de.x[0]]))
    np.savetxt('output_path.txt', alfa_v, fmt='%.3f', delimiter=',', header='Best PATH', comments='')
    print(alfa_v)
    spline, derivative, angle, angle_derivative, angle_sec_derivative = \
            path_info(left, right, alfa_v,N_angle)
    spline_points = spline(np.linspace(0,1,num = 1000))
    R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
    t1_SOCP_b,decision_variables_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)
    
    

    #Animates initial guess vs optimized solution
    animation_complete(spline,right,left,alfa_v,spline_points,decision_variables_SOCP_b,\
               t1_SOCP_b,n_discretization,m,mu,n_wheels)