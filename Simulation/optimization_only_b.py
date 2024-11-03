import numpy as np
import scipy as scp


# Create the vectors F_it for both the objective and constraints
# Input R_t (2d matrix with n_discretization vector), M_t (2d array with n_discretization vector M_t)
# C_tM_t (2d array with n_discretization vector C_t), n_discretization
# Output F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t)

def create_F_it(R_t,M_t,C_t,n_discretization):
    
    #create the n_dicretization vectors F_1t and F_2t
    F_1t, F_2t = np.zeros(np.shape(M_t)), np.zeros(np.shape(M_t))
    
    #loop to build the vectors
    for i in range(n_discretization-1):
        F_1t[i] = np.linalg.inv(R_t[i])@(M_t[i]/(2*1/(n_discretization-1))+C_t[i]/2)
        F_2t[i] = np.linalg.inv(R_t[i])@(-M_t[i]/(2*1/(n_discretization-1))+C_t[i]/2)
    return F_1t, F_2t




#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio vector A_t), 
#F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t), number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t, F_1t,F_2t,n_discretization,expansion_factor):
    def objective_function(b):
        decision_variables = b/expansion_factor
        cost=0
        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+2*xsi/(decision_variables[i+1]**0.5+decision_variables[i]**0.5)+(1-xsi)*np.transpose(decision_variables[i+1]*F_1t[i]+decision_variables[i]*F_2t[i])@A_t[i]
        return cost
    return objective_function


"""
#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,F_1t, F_2t,n_discretization):
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    def constraint1(decision_variables):
        remainder = np.zeros(2*(n_discretization-1))
        for i in range(n_discretization-1):
            val = decision_variables[i+1]*F_1t[i]+decision_variables[i]*F_2t[i]
            #first Cartesian coordinate 
            remainder[i]=decision_variables[u1+i]-val[0]
            #second Cartesian coordiinate
            remainder[n_discretization-1+i]=decision_variables[u2+i] - val[1]
        return remainder
    return constraint1
"""

#creates bounds to b, forcing it to be positive
def create_b_bounds(n_discretization):
    lb=[]
    ub=[]
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    #lb.extend([-np.inf]*2*(n_discretization-1))
    ub.extend([np.inf]*(n_discretization))#ub.extend([np.inf]*(3*n_discretization-2))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds

#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint(mu,mass, F_1t, F_2t, n_discretization,expansion_factor):
    def constraint(b):
        decision_variables = b/expansion_factor
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = mu*mass*9.81-np.linalg.norm((decision_variables[i+1]*F_1t[i]+decision_variables[i]*F_2t[i]))
        return remainder
    return constraint

"""
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
"""

#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_only_b(R_t,M_t,C_t,A_t,n_discretization,xsi,display):
    
    expansion_factor = 1E3
    
    
    #Creating force matrices F_1t and F_2t
    F_1t, F_2t = create_F_it(R_t,M_t,C_t,n_discretization)
    
    #creating objective and constraints
    objective_function = create_objective(xsi,A_t, F_1t, F_2t,n_discretization,expansion_factor)
    #constraint1 = create_constraint1(R_t,M_t,C_t,F_1t,F_2t,n_discretization)
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint =create_constraint(mu,mass,F_1t, F_2t, n_discretization,expansion_factor)
    bounds = create_b_bounds(n_discretization)
    
    
    cons = [
    #{'type': 'eq', 'fun': constraint1},  # Equality constraint 1
    {'type': 'ineq', 'fun': constraint}  # Inequality friction circle
        ]
    
    #optimizer options
    options = {
    'disp': display,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
    'ftol': 1e-8,      # Tolerance on function value changes
        }   
    
    # def callback_func(xk):
    #     callback_func.iteration += 1
    #     #print(f"Iteration {callback_func.iteration}")
    # callback_func.iteration = 0
    

    #creates initial guess inside the friction circle 
    x0 = np.ones(n_discretization)*expansion_factor
    
    # while not all sections forces inside the friction circle, reduce spline velocity in half
    while not ((constraint(x0)>= -1E-6).all()):
        x0=x0/2
    
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    if display:
        print("Test friction circle", (constraint(result.x)>= -1E-6).all())
        print("Test friction circle initial guess", (constraint(x0)>= -1E-6).all())
    result.x = result.x/expansion_factor
    return  result, x0




















"""
# Create the vectors F_it for both the objective and constraints
# Input R_t (2d matrix with n_discretization vector), M_t (2d array with n_discretization vector M_t)
# C_tM_t (2d array with n_discretization vector C_t), n_discretization
# Output F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t)
def create_F_it(R_t,M_t,C_t,n_discretization):
    
    #create the n_dicretization vectors F_1t and F_2t
    F_1t, F_2t = np.zeros(np.shape(M_t)), np.zeros(np.shape(M_t))
    
    #loop to build the vectors
    for i in range(n_discretization-1):
        F_1t[i] = np.linalg.inv(R_t[i])@(C_t[i]/2+M_t[i]/(2*1/(n_discretization-1)))
        F_2t[i] = np.linalg.inv(R_t[i])@(C_t[i]/2-M_t[i]/(2*1/(n_discretization-1)))
    return F_1t, F_2t


#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio vector A_t), 
#F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t), number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t, F_1t, F_2t,n_discretization):
    def objective_function(decision_variables):
        cost=0
        #sum of the objective over the path
        for i in range(n_discretization-1):
            cost = cost+(2*xsi)/(decision_variables[i+1]**0.5+decision_variables[i]**0.5)
            +(1-xsi)*np.transpose(decision_variables[i+1]*F_1t[i]+decision_variables[i]*F_2t[i])@A_t[i]
        return cost
    return objective_function


#creates bounds to b, forcing it to be positive
def create_b_bounds(n_discretization):
    lb=[]
    ub=[]
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    ub.extend([np.inf]*n_discretization)
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds

#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint(mu,mass,F_1t,F_2t,n_discretization):
    def constraint(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = mu*mass*9.81-np.linalg.norm(decision_variables[i+1]*F_1t[i]+decision_variables[i]*F_2t[i])
        print("Test friction circle", (remainder>= -1E-6).all())
        return remainder
    return constraint


#defines innequality constraint jacobian (friction circle)
#Input force equivalence vectors
#Output jacobian matrix
def create_constraint_jac(F_1t,F_2t):
    def constraint_jac(b):
        jac = np.zeros((len(b)-1,len(b)))
        #loop over the n-1 midpoints restrictions
        for i in range(len(b)-1):
            jac[i][i] = -(2*b[i]*np.transpose(F_2t[i])@F_2t[i]+b[i+1]*(np.transpose(F_1t[i])@F_2t[i]+np.transpose(F_2t[i])@F_1t[i]))
            jac[i][i+1] = -(2*b[i+1]*np.transpose(F_1t[i])@F_1t[i]+b[i]*(np.transpose(F_1t[i])@F_2t[i]+np.transpose(F_2t[i])@F_1t[i]))
        return jac
    
    return constraint_jac


def calc_constraint_jac(constraint, x, eps=1e-8):
    return scp.optimize.approx_fprime(x, constraint, eps)


#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_only_b(R_t,M_t,C_t,A_t,n_discretization,xsi):
    #building innitial guess
    x0 = np.ones(n_discretization)*0.0001
    
    #Creating force matrices F_1t and F_2t
    F_1t, F_2t = create_F_it(R_t,M_t,C_t,n_discretization)
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    #creating objective and constraints
    objective_function = create_objective(xsi, A_t, F_1t, F_2t,n_discretization)
    constraint = create_constraint(mu,mass,F_1t,F_2t,n_discretization)
    constraint_jac = create_constraint_jac(F_1t, F_2t)
    
    
    bounds = create_b_bounds(n_discretization)
    
    cons = [
    {'type': 'ineq', 'fun': constraint}  # Inequality friction circle , 'jac': constraint_jac 
        ]
    
    #optimizer options
    options = {
    'disp': True,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
    'ftol': 1e-10,     # Tolerance on function value changes
    'eps': 1e-12
        }   
    
    
    def callback_func(xk):
        callback_func.iteration += 1
        print(f"Iteration {callback_func.iteration}")

    callback_func.iteration = 0
    
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', constraints=cons,bounds=bounds,options=options, callback = callback_func)
    print("Test friction circle", (constraint(result.x)>= -1E-6).all())
    print("Test friction circle initial condition", (constraint(x0)>= -1E-6).all())
    print("Test b positive", (result.x>= 0).all())
    #print("Jacobiano da solução:", calc_constraint_jac(constraint,result.jac))
    #print("Jacobiano da solução calc:", constraint_jac(result.x))
    #if not (constraint(result.x)>= -1E-6).all():
    #    print("min", min(constraint(result.x)), constraint(result.x)) 
    #    print(constraint(x0))
     #   print(x0)
    return  result, x0
"""