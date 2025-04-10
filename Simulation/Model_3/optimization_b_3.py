import numpy as np
import scipy as scp
import matplotlib as plt





# Create the vectors F_it for both the objective and constraints
# Input R_t (2d matrix with n_discretization vector), M_t (2d array 
# with n_discretization vector M_t)
# C_tM_t (2d array with n_discretization vector C_t), n_discretization
# Output F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t)

def create_F_it(R_t,M_t,C_t,n_discretization):
    
    #create the n_dicretization vectors F_1t and F_2t
    F_1t, F_2t = np.zeros((n_discretization-1,8)),\
        np.zeros((n_discretization-1,8))
    
    #loop to build the vectors
    for i in range(n_discretization-1):
        F_1t[i] = np.linalg.pinv(R_t[i])@(M_t[i]/(2*1/(n_discretization-1))+
                                         C_t[i]/2)
        F_2t[i] = np.linalg.pinv(R_t[i])@(-M_t[i]/(2*1/(n_discretization-1))+
                                         C_t[i]/2)
    return F_1t, F_2t
















#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio 
# vector A_t), F_1t and F_2t (2d array with n_discretization vector F_1t and 
#vector F_2t), number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,F_1t,F_2t,T0,E0,n_discretization,expansion_factor):
    
    def objective_function(b):
        
        #expand space to facilitate the solver
        decision_variables = b/expansion_factor
        cost=0
        
        #sum over the path 
        for i in range(n_discretization-1):
            
            cost = cost + 2*xsi/((decision_variables[i+1]**0.5+decision_variables[i]
                    **0.5))+(1-xsi)*np.transpose(decision_variables[i+1]*F_1t[i]+
                                          decision_variables[i]*F_2t[i])@A_t[i]
                
        return cost
    
    return objective_function












#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,F_1t,F_2t,T0,E0,n_discretization,n_wheels,expansion_factor):
    def Grad(t):
        
        #expand space to facilitate the solver
        x= t/expansion_factor
        
        f = np.zeros(n_discretization)

        for i in range(n_discretization-1):
            f[i] += 2*xsi*(-1/(2*np.sqrt(x[i])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0+\
                (1-xsi)/E0*np.transpose(F_2t[i])@A_t[i]
            f[i+1] += 2*xsi*(-1/(2*np.sqrt(x[i+1])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0+\
                (1-xsi)/E0*np.transpose(F_1t[i])@A_t[i]
        return f/expansion_factor
    return Grad











#creates bounds to b 
def create_b_bounds(n_discretization):
    
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    ub.extend([np.inf]*(n_discretization))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds









#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint(mu,mass, F_1t, F_2t,n_discretization,n_wheels,expansion_factor):
    

    def constraint(b):
        decision_variables = b/expansion_factor
        remainder = np.zeros(n_wheels*(n_discretization-1))
        
        #in each section
        for i in range(n_discretization-1):
            #for every wheel
            for j in range(n_wheels):
                remainder[i+j*(n_discretization-1)] =mu*mass*9.81/n_wheels-\
                    np.linalg.norm((decision_variables[i+1]
                *F_1t[i][2*j:2*j+2]+decision_variables[i]*F_2t[i][2*j:2*j+2]))
        return remainder
    return constraint














#creates friction circle constraints gradient
def create_constraint_jac(F_1t, F_2t,n_discretization,n_wheels,expansion_factor):
    def constraint_jac(t):
         #expand space to facilitate the solver
        x = t/expansion_factor
        B1=np.zeros((n_wheels*(n_discretization-1),n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            for j in range(n_wheels):
                v = x[i+1]*F_1t[i][2*j:2*j+2]+x[i]*F_2t[i][2*j:2*j+2]
                norm_v = np.linalg.norm(v)
                if norm_v<1e-10:
                    norm_v=1e-10
                B1[i+j*(n_discretization-1),i] = -np.dot(F_2t[i][2*j:2*j+2],v)/norm_v
                B1[i+j*(n_discretization-1),i+1] = -np.dot(F_1t[i][2*j:2*j+2],v)/norm_v
        return B1/expansion_factor
    return constraint_jac










#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass 
# and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors 
# A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_b_3(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):
    if n_wheels != 4:
        print("Wrong optimization model. This one is specific for model3 (4 wheels)")
        SystemExit
    
    expansion_factor = 1E4
    T0=1
    E0=1
    
    #Creating force matrices F_1t and F_2t
    F_1t, F_2t = create_F_it(R_t,M_t,C_t,n_discretization)
    
    #creating objective and constraints
    objective_function = create_objective(xsi,A_t, F_1t, F_2t,\
        T0,E0,n_discretization,expansion_factor)
    grad = create_gradient_objective(xsi,A_t,F_1t,F_2t,T0,E0,n_discretization,n_wheels,expansion_factor)
    
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint =create_constraint(mu,mass,F_1t, F_2t, n_discretization,
                                  n_wheels,expansion_factor)
    constraint_jac =create_constraint_jac(F_1t, F_2t, n_discretization,
                                  n_wheels,expansion_factor)
    bounds = create_b_bounds(n_discretization)
    
    cons = [
    {'type': 'ineq', 'fun': constraint,'jac':constraint_jac}  # Inequality friction circle
        ]
    
    #optimizer options
    options = {
    'disp': display,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
    'ftol': 1e-8      # Tolerance on function value changes
        }   
    
    # def callback_func(xk):
    #     callback_func.iteration += 1
    #     #print(f"Iteration {callback_func.iteration}")
    # callback_func.iteration = 0
    
    #creates initial guess inside the friction circle 
    x0 = np.ones(n_discretization)*expansion_factor
    
    # while not all sections forces inside the friction circle, reduce 
    # spline velocity in half
    while not ((constraint(x0)>= -1E-6).all()):
        x0=x0/2

 
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP'
                        ,jac=grad,constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    
    
    if display:
        print("Test friction circle ", (constraint(result.x)>= -1E-6).all())
        print("Test friction circle initial guess ", (constraint(x0)>= -1E-6).all())
        
    result.x = result.x/expansion_factor
    return  result, x0