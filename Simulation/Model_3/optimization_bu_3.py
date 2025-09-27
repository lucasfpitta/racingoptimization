import numpy as np
import scipy as scp
import matplotlib as plt


#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization, number of wheels
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization
    
    def objective_function(t):
        #expand space to facilitate the solver
        decision_variables = np.zeros(len(t))
        decision_variables[0:n_discretization] = t[0:n_discretization]/expansion_factor_ab
        decision_variables[n_discretization:len(t)] = t[n_discretization:len(t)]/expansion_factor_u
        
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











#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization

    
    def Grad(t):
        
        #expand space to facilitate the solver
        x = np.zeros(len(t))
        x[0:2*n_discretization-1] = t[0:2*n_discretization-1]/expansion_factor_ab
        x[2*n_discretization-1:len(t)] = t[2*n_discretization-1:len(t)]/expansion_factor_u
        
        f = np.zeros((2*n_wheels)*(n_discretization-1)+n_discretization)
        
        #for each path section
        for i in range(n_discretization-1):
            for j in range(2*n_wheels):
                f[u+i+j*(n_discretization-1)]=(1-xsi)*A_t[i][j]/E0/expansion_factor_u

        for i in range(n_discretization-1):
            f[i] += 2*xsi*(-1/(2*np.sqrt(x[i])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0/expansion_factor_ab
            f[i+1] += 2*xsi*(-1/(2*np.sqrt(x[i+1])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0/expansion_factor_ab
        return f
    return Grad













#defines first equality constraint (dynamics)
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), number 
#of discretization
#Output 1d vector remainder, which goes to zero when the equality holds
def create_constraint1(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization
    
    def constraint1(t):
        
        #expand space to facilitate the solver
        decision_variables = np.zeros(len(t))
        decision_variables[0:n_discretization] = t[0:n_discretization]/expansion_factor_ab
        decision_variables[n_discretization:len(t)] = t[n_discretization:len(t)]/expansion_factor_u
        
        #Remainder array tells if the constraint is respected
        remainder = np.zeros(3*(n_discretization-1))
        
        for i in range(n_discretization-1):
            
            #first Cartesian coordinate 
            remainder[i] = -M_t[i][0]*((decision_variables[i+1]\
                        -decision_variables[i])/(2*1/(n_discretization-1)))-(
                    C_t[i][0])*(decision_variables[i]+decision_variables[i+1])/2
            
            #second Cartesian coordinate
            remainder[n_discretization-1+i] = -M_t[i][1]*((decision_variables[i+1]
                        -decision_variables[i])/(2*1/(n_discretization-1)))-(
                    C_t[i][1])*(decision_variables[i]+decision_variables[i+1])/2
            
            #angle coordinate
            remainder[2*(n_discretization-1)+i] = -M_t[i][2]*\
                (decision_variables[i+1]-decision_variables[i])/(2*1/(n_discretization-1))\
                    -(C_t[i][2])*(decision_variables[i]\
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













#defines equality constraint matrix F
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal 
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), 
# number of discretization
#Output constraint Matrix F
def create_constraint1_jac(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization
    
    F = np.zeros((3*(n_discretization-1),(2*n_wheels)*(n_discretization-1)+n_discretization))
        
        
        #iterate over each section to have the dynamics constraint on the 
        #two cartesian coordinates 
        #and the differential constraint
    for i in range(n_discretization-1):
            
            #first Cartesian coordinate dynamic constraint
        F[i,i]=(-C_t[i][0]/2+M_t[i][0]/(2*1/(n_discretization-1)))/expansion_factor_ab
        F[i,i+1]=(-C_t[i][0]/2-M_t[i][0]/(2*1/(n_discretization-1)))/expansion_factor_ab
        for j in range(2*n_wheels):
            F[i,u+i+j*(n_discretization-1)]=R_t[i][0][j]/expansion_factor_u  
            
            
            #second Cartesian coordinate dynamic constraint
        F[n_discretization-1+i,i]=(-C_t[i][1]/2+M_t[i][1]/(2*1/(n_discretization-1)))/expansion_factor_ab
        F[n_discretization-1+i,i+1]=(-C_t[i][1]/2-M_t[i][1]/(2*1/(n_discretization-1)))/expansion_factor_ab
        for j in range(2*n_wheels):
            F[n_discretization-1+i,u+i+j*(n_discretization-1)]=R_t[i][1][j]/expansion_factor_u  
                
                
                
            #Third Cartesian coordinate dynamic constraint
        F[2*(n_discretization-1)+i,i]=(-C_t[i][2]/2+M_t[i][2]/(2*1/(n_discretization-1)))/expansion_factor_ab
        F[2*(n_discretization-1)+i,i+1]=(-C_t[i][2]/2-M_t[i][2]/(2*1/(n_discretization-1)))/expansion_factor_ab
        for j in range(2*n_wheels):
            F[2*(n_discretization-1)+i,u+i+j*(n_discretization-1)]=R_t[i][2][j]/expansion_factor_u    

            
            
    def constraint1_jac(j):
        return F
    return constraint1_jac




















#creates bounds to b 
def create_b_bounds(n_discretization,n_wheels):
    
    length_wheels = 2*n_wheels*(n_discretization-1)
    
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*(length_wheels))
    ub.extend([np.inf]*(n_discretization+length_wheels))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds









#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_constraint2(mu,mass,n_discretization,n_wheels,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization

    def constraint2(t):
        decision_variables = t/expansion_factor_u
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
    return constraint2











#creates friction circle constraints gradient
def create_constraint2_jac(n_discretization,n_wheels,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization
    
    def constraint2_jac(t):
        #expand space to facilitate the solver
        x = t/expansion_factor_u
        B1=np.zeros((n_wheels*(n_discretization-1),\
            (2*n_wheels)*(n_discretization-1)+n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            for j in range(n_wheels):
                norm = np.sqrt(x[u+i+2*j*(n_discretization-1)]**2+x[u+i+(2*j+1)*(n_discretization-1)]**2)
                if norm > 1e-10:
                    B1[i+j*(n_discretization-1),u+i+2*j*(n_discretization-1)] = -x[u+i+2*j*(n_discretization-1)]/norm
                    B1[i+j*(n_discretization-1),u+i+(2*j+1)*(n_discretization-1)] = -x[u+i+(2*j+1)*(n_discretization-1)]/norm
                else: 
                    B1[i+j*(n_discretization-1),u+i+2*j*(n_discretization-1)] = -x[u+i+2*j*(n_discretization-1)]/1e-10
                    B1[i+j*(n_discretization-1),u+i+(2*j+1)*(n_discretization-1)] = -x[u+i+(2*j+1)*(n_discretization-1)]/1e-10
        return B1/expansion_factor_u
    return constraint2_jac





#Helps building innitial guess of constant b and normalization T0, E0
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal 
#C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #creates innitial guess
    x0 = (np.ones(n_discretization))*b0*expansion_factor_ab
    x0 = np.append(x0,np.zeros((2*n_wheels)*(n_discretization-1)))
    
    #flattened vector coordinates first 3 forces
    u1=n_discretization
    u2=u1+n_discretization-1
    u3=u2+n_discretization-1
    
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        u = b0*np.linalg.pinv(R_t[i])@C_t[i]*expansion_factor_u
            
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
def optimization_bu_3(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):
    if n_wheels != 3:
        print("Wrong optimization model. This one is specific for model3 (4 wheels)")
        SystemExit
    
    
    expansion_factor_ab = 1
    expansion_factor_u = 1
    
    
    #creating constraints
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    constraint1_jac = create_constraint1_jac(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    constraint2 =create_constraint2(mu,mass,n_discretization,n_wheels,expansion_factor_u)
    constraint2_jac =create_constraint2_jac(n_discretization,n_wheels,expansion_factor_u)
    bounds = create_b_bounds(n_discretization,n_wheels)
    
    cons = [
    {'type': 'eq', 'fun': constraint1,'jac': constraint1_jac},  # Equality constraint 1
    {'type': 'ineq', 'fun': constraint2,'jac': constraint2_jac}  # Inequality friction circle
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
    x0 , T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    while not ((constraint2(x0)>= -1E-6).all()):
        b0=b0/2
        x0, T0, E0 =  build_x0(b0,R_t,M_t,C_t,A_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
     #creating constraints

    E0=1
    T0=1
    objective_function = create_objective(xsi, A_t,abs(T0),abs(E0),\
        n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    grad = create_gradient_objective(xsi, A_t,abs(T0),abs(E0),\
        n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
 
    #optimization    
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP'
                        ,jac = grad,constraints=cons,bounds=bounds,options=options)#, callback = callback_func
    
    
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (constraint2(result.x)>= -1E-6).all())
        print("Test friction circle initial guess ", (constraint2(x0)>= -1E-6).all())
    result.x/=expansion_factor_ab    
    return  result, x0