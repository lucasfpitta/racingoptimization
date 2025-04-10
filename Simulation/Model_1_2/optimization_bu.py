import numpy as np
import scipy as scp




#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with
# n_discretizatio vector A_t),number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t,n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    
    
    def objective_function(decision_variables):
        cost=0
        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+2*xsi/(decision_variables[i+1]**0.5+decision_variables[i]
                    **0.5)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+
                                    decision_variables[u2+i]*A_t[i][1])
                    
        return cost
    return objective_function




#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    b=0
    u1=n_discretization
    u2=u1+n_discretization-1
    
    def Grad(x):

        f = np.zeros(2*(n_discretization-1)+n_discretization)
        
        #for each path section
        for i in range(n_discretization-1):
            f[u1+i]=(1-xsi)*A_t[i][0]/E0
            f[u2+i]=(1-xsi)*A_t[i][1]/E0

        for i in range(n_discretization-1):
            f[b+i] += 2*xsi*(-1/(2*np.sqrt(x[b+i])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0
            f[b+i+1] += 2*xsi*(-1/(2*np.sqrt(x[b+i+1])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0
        return f
    return Grad







#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), 
#number of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    
    
    def constraint1(decision_variables):
        
        remainder = np.zeros(2*(n_discretization-1))
        
        for i in range(n_discretization-1):
            #first Cartesian coordinate 
            remainder[i]=(R_t[i][0][0]*decision_variables[u1+i]+R_t[i][0][1]*
                        decision_variables[u2+i])-M_t[i][0]*((decision_variables[i+1]
                        -decision_variables[i])/(2*1/(n_discretization-1)))-(
                        C_t[i][0])*(decision_variables[i]+decision_variables[i+1])/2
                
            #second Cartesian coordiinate
            remainder[n_discretization-1+i]=(R_t[i][1][0]*decision_variables[u1+i]
                        +R_t[i][1][1]*decision_variables[u2+i])-M_t[i][1]*((
                        decision_variables[i+1]-decision_variables[i])/(2*1/(
                            n_discretization-1)))-C_t[i][1]*(decision_variables[i]+
                                                        decision_variables[i+1])/2
                            
        return remainder
    return constraint1















#defines first equality constraint (dynamics) jacobian
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number 
#of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1_jac(R_t,M_t,C_t,n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    
    #Remainder array tells if the constraint is respected
    F = np.zeros((2*(n_discretization-1),n_discretization+2*
                  (n_discretization-1)))

        
    for i in range(n_discretization-1):
            
        #first Cartesian coordinate dynamic constraint
        F[i,i]=-C_t[i][0]/2+M_t[i][0]/(2*1/(n_discretization-1))
        F[i,i+1]=-C_t[i][0]/2-M_t[i][0]/(2*1/(n_discretization-1))
        F[i,u1+i]=R_t[i][0][0]
        F[i,u2+i]=R_t[i][0][1]
            
            
            #second Cartesian coordinate dynamic constraint
        F[n_discretization-1+i,i]=-C_t[i][1]/2+M_t[i][1]/(2*1/(n_discretization-1))
        F[n_discretization-1+i,i+1]=-C_t[i][1]/2-M_t[i][1]/(2*1/(n_discretization-1))
        F[n_discretization-1+i,u1+i]=R_t[i][1][0]
        F[n_discretization-1+i,u2+i]=R_t[i][1][1]
        
    def constraint1_jac(decision_variables):         
        return F
    return constraint1_jac













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
            remainder[i] = -(decision_variables[u1+i]**2+decision_variables[u2+i]
                             **2)**0.5+mu*mass*9.81
        return remainder
    return constraint2





#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint2_jac(n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization
    u2=n_discretization+n_discretization-1
    
    
    def constraint2_jac(x):
        B1=np.zeros((n_discretization-1,n_discretization+2*(n_discretization-1)))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            norm = np.sqrt(x[u1 + i]**2 + x[u2 + i]**2)
            B1[i,u1+i] = -x[u1+i]/norm#/(x[u1+i]**2+x[u2+i]**2)**0.5
            B1[i,u2+i] = -x[u2+i]/norm#/(x[u1+i]**2+x[u2+i]**2)**0.5
        return B1
    return constraint2_jac










#Helps building innitial guess of constant b
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal  C_t 
#2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0, R_t,C_t,n_discretization):
    
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
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
# Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, M_t 
#and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_bu(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):
    
    E0=1
    T0=1
    
    #creating objective and constraints
    objective_function = create_objective(xsi, A_t,n_discretization)
    grad = create_gradient_objective(xsi,A_t,abs(T0),abs(E0),n_discretization)
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization)
    constraint1_jac = create_constraint1_jac(R_t,M_t,C_t,n_discretization)
    
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    
    constraint2 =create_constraint2(mu,mass,n_discretization)
    constraint2_jac =create_constraint2_jac(n_discretization)
    bounds = create_b_bounds(n_discretization)
    
    
    cons = [
    {'type': 'eq', 'fun': constraint1,'jac':constraint1_jac},  # Equality constraint 1
    {'type': 'ineq', 'fun': constraint2,'jac':constraint2_jac}  # Inequality friction circle
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
    
    
    b0 = 1
    #building innitial guess
    x0 =  build_x0(b0,R_t,C_t,n_discretization)
    while not ((constraint2(x0)>= -1E-6).all()):
        b0=b0/2
        x0 =  build_x0(b0,R_t,C_t,n_discretization)
    
    
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', 
                        jac=grad,constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    if display:
        print("Test friction circle", (constraint2(result.x)>= -1E-6).all())
        print("Test friction circle initial guess", (constraint2(x0)>= -1E-6).all())
    return  result, x0