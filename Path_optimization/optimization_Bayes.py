import numpy as np
from Map_processing.mapreader import get_outline
from splines.splines import path_info
from Physics.model2 import model2
from Simulation.optimization_main import init_optimization_SOCP_b
import matplotlib.pyplot as plt
from Visualization.plots import animation_complete
from scipy.optimize import differential_evolution
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from Physics.translate import translate_velocity



from skopt import gp_minimize
from skopt.space import Real




class LivePlotter:
    def __init__(self):
        self.iteration_data = []
        self.cost_data = []
        self.best_cost = []

        
        # Set up the plot
        plt.ion()  # Interactive mode
        self.fig, (self.ax1, self.ax2) = plt.subplots(ncols=2, figsize=(12, 6))
        self.line, = self.ax1.plot([], [], label="Cost", color='blue', linewidth=2)
        self.line2, = self.ax1.plot([], [], label="Best Cost", color='red', linewidth=2)
        self.ax1.set_title("Optimization Convergence")
        self.ax1.set_xlabel("Iteration")
        self.ax1.set_ylabel("Cost")
        self.ax1.legend()
        self.ax1.grid()

        self.alfas_data = []

        self.ax2.set_title("Best Trajectory")
        self.ax2.grid()
        self.ax2.autoscale()  # Automatically scale the plot to the data
        self.ax2.set_xlabel("X")
        self.ax2.set_ylabel("Y")
       
        self.norm = Normalize(vmin=0, vmax=1)
        self.cmap = plt.cm.RdYlGn_r


        # Add colorbar (only once)
        self.colorbar = plt.colorbar(
            plt.cm.ScalarMappable(norm=self.norm, cmap=self.cmap),
            ax=self.ax2,
            orientation="vertical"
        )

        self.colorbar.set_label("Velocity (m/s)")
        plt.tight_layout()  # Adjust layout to prevent overlap


    def update(self, iteration, cost, alfas, b,spline,derivative,n_discretization,left,right):
        # Update the data
        self.iteration_data.append(iteration)
        self.cost_data.append(cost)
        
        if len(self.best_cost)==0 or cost<self.best_cost[-1]:
            self.best_cost.append(cost)
            self.alfas_data.append(np.concatenate((alfas,[alfas[0]])))

            spline_points = spline(np.linspace(0,1,num = n_discretization))
            velocity = translate_velocity(derivative,b,n_discretization)
            points = np.array([spline_points[0], spline_points[1]]).T.reshape(-1, 1, 2) 
            segments = np.concatenate([points[:-1], points[1:]], axis=1)  # Create line segments

            self.norm = Normalize(vmin=np.min(velocity), vmax=np.max(velocity))


            if len(velocity) != len(spline_points[0]) - 1:
                print(f"Velocity length mismatch: {len(velocity)} vs {len(spline_points[0]) - 1}")
                return  # Abort if there's a mismatch

            for collection in list(self.ax2.collections):
                collection.remove()
            
            

            lc = LineCollection(segments, cmap=self.cmap, norm=self.norm)
            lc.set_array(velocity)  # Attach velocity data
            lc.set_linewidth(2)  # Set line width
            self.ax2.add_collection(lc)


            # Update the colorbar
            self.colorbar.mappable.set_array(velocity)  # Set the data for the colorbar
            self.colorbar.mappable.set_norm(self.norm)  # Update the norm
            self.colorbar.update_ticks()  # Update the ticks on the colorbar


            for line in self.ax2.lines:
                line.remove()
            #self.ax2.plot(spline_points[0], spline_points[1], color='blue')  # Simple plot
            self.ax2.plot(left[0], left[1], color='black',linewidth=1, alpha=0.5)  # Simple plot
            self.ax2.plot(right[0], right[1], color='black',linewidth=1, alpha=0.5)  # Simple plot


            # Adjust axis limits
            self.ax2.relim()  # Recalculate limits
            self.ax2.autoscale_view()  # Rescale view
            self.fig.canvas.draw_idle()  # Ensure drawing is updated
    

        else:
            self.best_cost.append(self.best_cost[-1])

        # Update the plot
        self.line.set_data(self.iteration_data, self.cost_data)
        self.line2.set_data(self.iteration_data, self.best_cost)
        self.ax1.relim()  # Recalculate limits
        self.ax1.autoscale_view()  # Rescale view
        
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()  # Force the update to reflect immediately
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
        t1_SOCP_b,b =init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)
        cost = t1_SOCP_b[-1]


        
        if plotter is not None:
            iteration = len(plotter.iteration_data)  # Determine current iteration
            plotter.update(iteration, cost, alfas, b,spline,derivative,n_discretization,left,right)
        
        if len(plotter.best_cost)!=0:
            print(f"Current time: {cost:.2f}, Best time: {plotter.best_cost[-1]:.2f}")
        
        return cost
    return path_optimization







# Function to generate bounded normal perturbations
def generate_bounded_population(best_solution, bounds, population_size, std_dev=0.025):
    population = []
    for _ in range(population_size):
        individual = []
        for i, (low, high) in enumerate(bounds):
            # Generate a normally distributed value and clip it to the bounds
            value = np.clip(np.random.normal(best_solution[i], std_dev), low, high)
            individual.append(value)
        population.append(individual)
    return np.array(population)















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
    initial_population = generate_bounded_population(result_b.x, new_bounds, population_size)


    

    # Run Differential Evolution
    result_de = differential_evolution(
        func=wrapped_objective,
        bounds=new_bounds,
        #init=initial_population,  # Use the valid population format
        strategy="best2bin",
        mutation=(0.5, 0.9),  # Mutation range for more diverse steps
        recombination=0.7,    # Default recombination
        maxiter=100,
        init=initial_population,  # Parse the custom population
        tol=1e-6,
        seed=42
    )

    # Print the results
    print("Best score Differential Evolution:", result_de.fun)
    print("Best solution Differential Evolution:", result_de.x)

    plotter.finalize()
    alfa_v = np.concatenate((result_b.x,[result_b.x[0]]))
    np.savetxt('output_path.txt', alfa_v, fmt='%.3f', delimiter=',', header='Best PATH', comments='')


    spline, derivative, angle, angle_derivative, angle_sec_derivative = \
            path_info(left, right, alfa_v,N_angle)
    spline_points = spline(np.linspace(0,1,num = n_discretization))
    R_t, M_t, C_t, A_t = model2(spline,angle,n_discretization,m,mu,pho_air,A0,Cx)
    t1_SOCP_b,decision_variables_SOCP_b=init_optimization_SOCP_b(R_t, M_t, C_t, A_t,n_discretization,xsi,n_wheels,display=False,plot=True)
    
    

    #Animates initial guess vs optimized solution
    animation_complete(spline,right,left,alfa_v,spline_points,decision_variables_SOCP_b,\
               t1_SOCP_b,n_discretization,m,mu,n_wheels)
    

    velocity = translate_velocity(derivative,decision_variables_SOCP_b,n_discretization)


    # Prepare segments for LineCollection
    points = np.array([spline_points[0], spline_points[1]]).T.reshape(-1, 1, 2)  # Reshape to (N, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)  # Create line segments

    # Normalize the velocity for colormap
    norm = Normalize(vmin=np.min(velocity), vmax=np.max(velocity))
    cmap = plt.cm.RdYlGn_r  # Choose a colormap (e.g., green to red)

    # Create the LineCollection
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(velocity)  # Attach velocity data
    lc.set_linewidth(2)  # Set line width

    # Plot
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.autoscale()  # Automatically scale the plot to the data
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Trajectory Colored by Velocity")

    # Add colorbar
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.RdYlGn), ax=ax, orientation='vertical')
    cb.set_label("Velocity (m/s)")

    plt.show()