import numpy as np
import scipy as scp


#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,T0,E0,n_discretization,n_wheels):
    
    #flattened vector coordinates
    u=n_discretization+n_discretization-1
    
    def objective_function(decision_variables):
        
        cost=0
        
        #sum over the path 
        for i in range(n_discretization-1):
            
            cost = cost + 2*xsi/((decision_variables[i+1]**0.5+decision_variables[i]
                    **0.5)*T0)
            
            for j in range(2*n_wheels):
                cost = cost+(1-xsi)*(decision_variables[\
                    u+i+j*(n_discretization-1)]*A_t[i][j])/E0
                
        return cost
    
    return objective_function










#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number 
#of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,n_discretization,n_wheels):
    
    #flattened vector coordinates
    a = n_discretization
    u=n_discretization+n_discretization-1
    
    def constraint1(decision_variables):
        
        #Remainder array tells if the constraint is respected
        remainder = np.zeros(3*(n_discretization-1))
        
        for i in range(n_discretization-1):
            
            #first Cartesian coordinate 
            remainder[i] = -M_t[i][0]*decision_variables[a+i]-(
                    C_t[i][0])*(decision_variables[i]+decision_variables[i+1])/2
            
            #second Cartesian coordinate
            remainder[n_discretization-1+i] = -M_t[i][1]*decision_variables[a+i]-(
                    C_t[i][1])*(decision_variables[i]+decision_variables[i+1])/2
            
            #angle coordinate
            remainder[2*(n_discretization-1)+i] = -M_t[i][2]*\
                decision_variables[a+i]-(C_t[i][2])*(decision_variables[i]\
                    +decision_variables[i+1])/2
            
            for j in range(2*n_wheels):

                #first Cartesian coordinate 
                remainder[i]=remainder[i]+(R_t[i][0][j]*decision_variables[\
                    u+i+j*(n_discretization-1)])
                
                #second Cartesian coordinate
                remainder[n_discretization-1+i]=remainder[n_discretization-1+i]+\
                    (R_t[i][1][j]*decision_variables[u+i+j*(n_discretization-1)])
                        
                #angle coordinate
                remainder[2*(n_discretization-1)+i]=remainder[2*(\
                    n_discretization-1)+i]+(R_t[i][2][j]*decision_variables[u+i+j\
                        *(n_discretization-1)])
                    
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
            remainder[i] = decision_variables[i+1]-decision_variables[i]-(2
                            )*decision_variables[a+i]*1/(n_discretization-1)
        return remainder
    return constraint2









#creates bounds to b 
def create_b_bounds(n_discretization,n_wheels):
    
    length_wheels = 2*n_wheels*(n_discretization-1)
    
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*(n_discretization-1+length_wheels))
    ub.extend([np.inf]*(2*n_discretization-1+length_wheels))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds









#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint3(mu,mass,n_discretization,n_wheels):
    
    #flattened vector coordinates
    u=n_discretization+n_discretization-1

    def constraint3(decision_variables):
        remainder = np.zeros(n_wheels*(n_discretization-1))
        
        #in each section
        for i in range(n_discretization-1):
            
            #for every wheel
            for j in range(n_wheels):
                remainder[i+j*(n_discretization-1)] = -(\
                    decision_variables[u+i+2*j*(n_discretization-1)]**2\
                        +decision_variables[u+i+(2*j+1)*(n_discretization-1)]
                             **2)**0.5+mu*mass*9.81/n_wheels
            
        return remainder
    return constraint3









#Helps building innitial guess of constant b and normalization T0, E0
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal 
#C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels):
    
    #creates innitial guess
    x0 = (np.ones(n_discretization))*b0
    x0 = np.append(x0,np.zeros((1+2*n_wheels)*(n_discretization-1)))
    
    #flattened vector coordinates first 3 forces
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    u3=u2+n_discretization-1
    
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        a = (x0[i+1]-x0[i])/(2*1/(n_discretization-1))
        u = a*np.linalg.pinv(R_t[i])@M_t[i]+(x0[i+1]+x0[i])/2*\
            np.linalg.pinv(R_t[i])@C_t[i]
            
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]
        x0[u3+i]=u[2]
        
    #Calculates the normalization factors
    T0=0
    E0=0
    # for i in range(n_discretization-1):
    #    T0 = T0+2/(x0[i+1]**0.5+x0[i]**0.5)
    #    E0 = E0+(x0[u1+i]*A_t[i][0]+x0[u2+i]*A_t[i][1])
    return x0, T0, E0












#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass 
# and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors 
# A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_abu_3(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):
    if n_wheels != 4:
        print("Wrong optimization model. This one is specific for model3 (4 wheels)")
        SystemExit
    
    #creating constraints
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization,n_wheels)
    constraint2=create_constraint2(n_discretization)
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint3 =create_constraint3(mu,mass,n_discretization,n_wheels)
    bounds = create_b_bounds(n_discretization,n_wheels)
    
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
    x0 , T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels)
    while not ((constraint3(x0)>= -1E-6).all()):
        b0=b0/2
        x0, T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels)
     #creating constraints

    E0=1
    T0=1
    objective_function = create_objective(xsi, A_t,abs(T0),abs(E0),\
        n_discretization,n_wheels)
 
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP'
                        , constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    
    
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (constraint3(result.x)>= -1E-6).all())
        print("Test friction circle initial guess ", (constraint3(x0)>= -1E-6).all())
        
    return  result, x0