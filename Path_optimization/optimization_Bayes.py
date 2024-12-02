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




class BestResultTracker:
    def __init__(self):
        self.best_value = 1e5  # Start with a high value for minimization
        self.best_params = np.ones(8)         # To store the best parameters

    def update(self, alfas, cost):
        """
        Updates the best result if the current cost is better.
        :param alfas: Current parameters being evaluated
        :param cost: Current objective function value (cost)
        """
        print(self.best_value)
        print(cost)
        print(self.best_value is None or cost < self.best_value)
        if self.best_value is None or cost < self.best_value:
            print(f"New Best value: {cost}")
            self.best_value = cost
            self.best_params = alfas

    def get_best(self):
        """
        Returns the best result so far.
        :return: (best_params, best_value)
        """
        return self.best_params, self.best_value











def init_path_optimization(right,left,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels):

    def path_optimization(alfas):
        alfa_v = np.concatenate((alfas,[alfas[0]]))
        spline, derivative, angle, angle_derivative, angle_sec_derivative = \
        path_info(left, right, alfa_v,N_angle)
        R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
        t1_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=False)
        cost = t1_SOCP_b[-1]
        print("current cost: ", cost)
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
    tracker = BestResultTracker()

    # Run Bayesian optimization
    result_b = gp_minimize(
        func=objective,  # Objective function to minimize
        dimensions=space,         # Search space
        acq_func="EI",            # Acquisition function: Expected Improvement
        n_calls=200,               # Number of evaluations
        random_state=42           # Seed for reproducibility
    )

    # Print the results
    print("Best score Bayes:", result_b.fun)
    print("Best solution Bayes:", result_b.x)

    new_bounds=[]
    for i in range(dimensions):
        new_bounds.append((max(result_b.x[i]-0.05,0),min(result_b.x[i]+0.05,1)))

     # Expand the Bayesian result to form a valid population
    population_size = 10  

    # Run Differential Evolution
    result_de = differential_evolution(
        func=objective,
        bounds=new_bounds,
        #init=initial_population,  # Use the valid population format
        strategy="best2bin",
        mutation=(0.5, 0.9),  # Mutation range for more diverse steps
        recombination=0.7,    # Default recombination
        maxiter=1000,
        popsize=population_size,  # Ensure consistency with init size
        tol=1e-6,
        seed=42
    )

    # Print the results
    print("Best score Differential Evolution:", result_de.fun)
    print("Best solution Differential Evolution:", result_de.x)


    alfa_v = np.concatenate((result_de.x,[result_de.x[0]]))
    print(alfa_v)
    spline, derivative, angle, angle_derivative, angle_sec_derivative = \
            path_info(left, right, alfa_v,N_angle)
    spline_points = spline(np.linspace(0,1,num = 1000))
    R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
    t1_SOCP_b,decision_variables_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)

    

    #Animates initial guess vs optimized solution
    animation_complete(spline,right,left,alfa_v,spline_points,decision_variables_SOCP_b,\
               t1_SOCP_b,n_discretization,m,mu,n_wheels)