import numpy as np
import scipy as scp
import matplotlib as plt

#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t,T0,E0,n_discretization):
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def objective_function(decision_variables):
        cost=0
        cost1=0
        cost2=0
        #sum over the path 
        for i in range(n_discretization-1):
            cost1=cost1+2/((decision_variables[i+1]**0.5+decision_variables[i]**0.5)*T0)
            cost2= cost2+(decision_variables[u1+i]*A_t[i][0]+decision_variables[u2+i]*A_t[i][1])/E0
            cost = cost+(2*xsi/((decision_variables[i+1]**0.5+decision_variables[i]**0.5)*T0)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+decision_variables[u2+i]*A_t[i][1])/E0)
        #print("Cost time: ", cost1, " Cost Energy: ", cost2)
        return cost
    
    return objective_function


#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,n_discretization):
    #flattened vector coordinates
    a = n_discretization
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def constraint1(decision_variables):
        remainder = np.zeros(2*(n_discretization-1))
        for i in range(n_discretization-1):
            #first Cartesian coordinate 
            remainder[i]=(R_t[i][0][0]*decision_variables[u1+i]+R_t[i][0][1]*decision_variables[u2+i])-M_t[i][0]*decision_variables[a+i]-C_t[i][0]*(decision_variables[i]+decision_variables[i+1])/2
            #second Cartesian coordiinate
            remainder[n_discretization-1+i]=(R_t[i][1][0]*decision_variables[u1+i]+R_t[i][1][1]*decision_variables[u2+i])-M_t[i][1]*decision_variables[a+i]-C_t[i][1]*(decision_variables[i]+decision_variables[i+1])/2
        return remainder
    return constraint1

#defines second equality constraint (differential)
#Input number of discretizations
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint2(n_discretization):
    #flattened vector coordinates
    a = n_discretization
    def constraint2(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = decision_variables[i+1]-decision_variables[i]-2*decision_variables[a+i]*1/(n_discretization-1)
        return remainder
    return constraint2

#creates bounds to b 
def create_b_bounds(n_discretization):
    lb=[]
    ub=[]
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-8]*n_discretization)
    lb.extend([-np.inf]*3*(n_discretization-1))
    ub.extend([1E-8])
    ub.extend([np.inf]*(4*n_discretization-4))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds

#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint3(mu,mass,n_discretization):
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def constraint3(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = -(decision_variables[u1+i]**2+decision_variables[u2+i]**2)**0.5+mu*mass*9.81
        return remainder
    return constraint3

#Helps building innitial guess of constant b and normalization T0, E0
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal  C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0,R_t,M_t,C_t, A_t,n_discretization):
    
    #creates innitial guess
    x0 = (np.ones(n_discretization)+np.random.uniform(low=-0.5, high=0.5, size=n_discretization))*b0
    x0 = np.append(x0,np.zeros(3*(n_discretization-1)))
    x0[0]=1E-8
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        a = (x0[i+1]-x0[i])/(2*1/(n_discretization-1))
        u = a*np.linalg.inv(R_t[i])@M_t[i]+(x0[i+1]+x0[i])/2*np.linalg.inv(R_t[i])@C_t[i]
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]
        
        
    #T0 = 20/(1/(n_discretization-1))
    #E0 = 1/2*85*(30/2.6)**2 
    T0=0
    E0=0
    for i in range(n_discretization-1):
       T0 = T0+2/(x0[i+1]**0.5+x0[i]**0.5)
       E0 = E0+(x0[u1+i]*A_t[i][0]+x0[u2+i]*A_t[i][1])
    return x0, T0, E0

#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization(R_t,M_t,C_t,A_t,n_discretization,xsi,display):
    
    #creating constraints
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization)
    constraint2=create_constraint2(n_discretization)
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint3 =create_constraint3(mu,mass,n_discretization)
    bounds = create_b_bounds(n_discretization)
    
    cons = [
    {'type': 'eq', 'fun': constraint1},  # Equality constraint 1
    {'type': 'eq', 'fun': constraint2},  # Equality constraint 2
    {'type': 'ineq', 'fun': constraint3}  # Inequality friction circle
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
    
    b0=1
    #building innitial guess
    x0 , T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization)
    while not ((constraint3(x0)>= -1E-6).all()):
        b0=b0/2
        x0, T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization)
    
     #creating constraints
    objective_function = create_objective(xsi, A_t,abs(T0),abs(E0),n_discretization)
 
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (constraint3(result.x)>= -1E-6).all())
        print("Test friction circle initial guess ", (constraint3(x0)>= -1E-6).all())
    return  result, x0