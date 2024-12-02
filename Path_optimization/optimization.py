import numpy as np
from Map_processing.mapreader import get_outline
from splines.splines import path_info
from Physics.model2 import model2
from Simulation.optimization_main import init_optimization_SOCP_b
import pyswarms as ps
from pyswarms.utils.functions import single_obj
from pyswarms.utils.plotters import plot_cost_history
import matplotlib.pyplot as plt
from Visualization.plots import animation_complete





class LivePlotter:
    def __init__(self):
        self.iteration_data = []
        self.best_cost_data = []
        
        # Set up the plot
        plt.ion()  # Interactive mode
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], label="Best Cost", color='blue', linewidth=2)
        self.ax.set_title("PSO Optimization Convergence")
        self.ax.set_xlabel("Iteration")
        self.ax.set_ylabel("Cost")
        self.ax.legend()
        self.ax.grid()

    def update(self, iteration, best_cost):
        # Update the data
        self.iteration_data.append(iteration)
        self.best_cost_data.append(best_cost)

        # Update the plot
        self.line.set_data(self.iteration_data, self.best_cost_data)
        self.ax.relim()  # Recalculate limits
        self.ax.autoscale_view()  # Rescale view
        plt.pause(0.01)  # Pause to refresh the plot

    def finalize(self):
        plt.ioff()
        plt.show()














def init_path_optimization(right,left,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels,n_particles):

    def path_optimization(alfas, plotter=None):
        #defines the ´path
        cost = np.zeros(n_particles)
        for i in range(n_particles):
            alfa_v = np.concatenate((alfas[i],[alfas[i][0]]))
            spline, derivative, angle, angle_derivative, angle_sec_derivative = \
            path_info(left, right, alfa_v,N_angle)
            R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
            t1_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=False)
            cost[i] = t1_SOCP_b[-1]
        print()
        print(np.min(cost))
        if plotter is not None:
            best_cost=np.min(cost)
            iteration = len(plotter.iteration_data)  # Determine current iteration
            plotter.update(iteration, best_cost)
        return cost
    return path_optimization


def trajectory_optimization(external,internal,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels):
    #defines the outline
    right, left,  = get_outline(external,internal)

    

    # Define problem dimensions and bounds
    dimensions = len(right[0])-1  # Number of control_points
    print(f"The number of control points is {dimensions}")
    lower_bound = 0
    upper_bound = 1
    bounds = (lower_bound * np.ones(dimensions), upper_bound * np.ones(dimensions))

    # Define PSO hyperparameters
    options = {
        'c1': 0.8,  # Cognitive parameter
        'c2': 1.5,  # Social parameter
        'w': 0.6    # Inertia weight
    }

    n_particl = 30
    # Initialize the PSO optimizer
    optimizer = ps.single.GlobalBestPSO(
        n_particles=n_particl,       # Number of particles
        dimensions=dimensions,
        options=options,
        bounds=bounds
    )

    objective = init_path_optimization(right,left,N_angle,n_discretization,m,mu,pho_air,A0,Cx,xsi,n_wheels,n_particl)
    
    plotter = LivePlotter()
    # Run optimization
    best_cost, best_position = optimizer.optimize(objective, iters=1500,plotter=plotter)

    plotter.finalize()
    # Print the results
    print(f"Best position: {best_position}")
    print(f"Best cost: {best_cost}")
    # Plot the convergence
    plot_cost_history(optimizer.cost_history)
    plt.title("Convergence Plot")
    plt.xlabel("Iterations")
    plt.ylabel("Cost")
    plt.show()

    alfa_v = np.concatenate((best_position,[best_position[0]]))
    print(alfa_v)
    spline, derivative, angle, angle_derivative, angle_sec_derivative = \
            path_info(left, right, alfa_v,N_angle)
    spline_points = spline(np.linspace(0,1,num = 1000))
    R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
    t1_SOCP_b,decision_variables_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)

    

    #Animates initial guess vs optimized solution
    animation_complete(spline,right,left,alfa_v,spline_points,decision_variables_SOCP_b,\
               t1_SOCP_b,n_discretization,m,mu,n_wheels)