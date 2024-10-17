import numpy as np
import scipy as scp

#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t,n_discretization):
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    def objective_function(decision_variables):
        cost=0
        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+2*xsi/(decision_variables[i+1]**0.5+decision_variables[i]**0.5)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+decision_variables[u2+i]*A_t[i][1])
        return cost
    return objective_function


#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,n_discretization):
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    def constraint1(decision_variables):
        remainder = np.zeros(2*(n_discretization-1))
        for i in range(n_discretization-1):
            #first Cartesian coordinate 
            remainder[i]=(R_t[i][0][0]*decision_variables[u1+i]+R_t[i][0][1]*decision_variables[u2+i])-M_t[i][0]*((decision_variables[i+1]-decision_variables[i])/(2*1/(n_discretization-1)))-C_t[i][0]*(decision_variables[i]+decision_variables[i+1])/2
            #second Cartesian coordiinate
            remainder[n_discretization-1+i]=(R_t[i][1][0]*decision_variables[u1+i]+R_t[i][1][1]*decision_variables[u2+i])-M_t[i][1]*((decision_variables[i+1]-decision_variables[i])/(2*1/(n_discretization-1)))-C_t[i][1]*(decision_variables[i]+decision_variables[i+1])/2
        return remainder
    return constraint1


#creates bounds to b 
def create_b_bounds(n_discretization):
    lb=[]
    ub=[]
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*2*(n_discretization-1))
    ub.extend([np.inf]*(3*n_discretization-2))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds

#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint2(mu,mass,n_discretization):
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    def constraint2(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = -(decision_variables[u1+i]**2+decision_variables[u2+i]**2)**0.5+mu*mass*9.81
        return remainder
    return constraint2

#Helps building innitial guess of constant b
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal  C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(R_t,C_t,n_discretization):
    #def constant b
    b0=0.0005
    #creates innitial guess
    x0 = np.ones(n_discretization)*b0
    x0 = np.append(x0,np.zeros(2*(n_discretization-1)))
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        u = b0*np.dot(np.linalg.inv(R_t[i]),C_t[i])
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]
    return x0

#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_bu(R_t,M_t,C_t,A_t,n_discretization,xsi):
    #building innitial guess
    x0 =  build_x0(R_t,C_t,n_discretization)#{"u": np.ones((n_discretization,2)), "a": np.ones(n_discretization-1), "b": np.ones(n_discretization-1) }
    
    #creating objective and constraints
    objective_function = create_objective(xsi, A_t,n_discretization)
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization)
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint2 =create_constraint2(mu,mass,n_discretization)
    bounds = create_b_bounds(n_discretization)
    
    cons = [
    {'type': 'eq', 'fun': constraint1},  # Equality constraint 1
    {'type': 'ineq', 'fun': constraint2}  # Inequality friction circle
        ]
    
    #optimizer options
    options = {
    'disp': True,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
    'ftol': 1e-8      # Tolerance on function value changes
        }   
    def callback_func(xk):
        callback_func.iteration += 1
        #print(f"Iteration {callback_func.iteration}")

    callback_func.iteration = 0
    
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', constraints=cons,bounds=bounds,options=options, callback = callback_func)
    print("Test friction circle", (constraint2(result.x)>= -1E-6).all())
    print("Test friction circle initial guess", (constraint2(x0)>= -1E-6).all())
    print("Test b positive", (result.x[0:n_discretization]>= 0).all())
    return  result, x0