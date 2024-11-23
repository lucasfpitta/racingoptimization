import numpy as np
import scipy as scp
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy import sparse 














#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    b=n_discretization-1
    u=b+n_discretization
    
    def objective_function(t):
        #expand space to facilitate the solver
        decision_variables = np.zeros(len(t))
        decision_variables[0:2*n_discretization-1] = t[0:2*n_discretization-1]/expansion_factor_ab
        decision_variables[2*n_discretization-1:len(t)] = t[2*n_discretization-1:len(t)]/expansion_factor_u
        
        cost=0
        
        #sum over the path 
        for i in range(n_discretization-1):
            
            cost = cost + 2*xsi/((decision_variables[b+i+1]**0.5+decision_variables[b+i]
                    **0.5)*T0)
            
            for j in range(2*n_wheels):
                cost = cost+(1-xsi)*(decision_variables[\
                    u+i+j*(n_discretization-1)]*A_t[i][j])/E0
                
        return cost
    
    return objective_function













#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    b=n_discretization-1
    u=b+n_discretization

    
    def Grad(t):
        
        #expand space to facilitate the solver
        x = np.zeros(len(t))
        x[0:2*n_discretization-1] = t[0:2*n_discretization-1]/expansion_factor_ab
        x[2*n_discretization-1:len(t)] = t[2*n_discretization-1:len(t)]/expansion_factor_u
        
        f = np.zeros((2*n_wheels+1)*(n_discretization-1)+n_discretization)
        
        #for each path section
        for i in range(n_discretization-1):
            for j in range(2*n_wheels):
                f[u+i+j*(n_discretization-1)]=(1-xsi)*A_t[i][j]/E0/expansion_factor_u

        for i in range(n_discretization-1):
            f[b+i] += 2*xsi*(-1/(2*np.sqrt(x[b+i])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0/expansion_factor_ab
            f[b+i+1] += 2*xsi*(-1/(2*np.sqrt(x[b+i+1])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0/expansion_factor_ab
        return f
    return Grad






#defines Hessian of the function at a linearization point
def create_Hessian_objective(xsi,T0,n_discretization,n_wheels,expansion_factor_ab):
    
    #flattened vector coordinates
    b=n_discretization-1
    
    def Hessian(t):
        #expand space to facilitate the solver
        x = t/expansion_factor_ab
    
        H = sparse.lil_matrix(((2*n_wheels+1)*(n_discretization-1)+n_discretization,\
            (2*n_wheels+1)*(n_discretization-1)+n_discretization))
        
        
        #for each path section
        for i in range(n_discretization-1):
            H[b+i,b+i]+=2*xsi*(x[b+i+1]**0.5+3*x[b+i]**0.5)/(4*x[b+i]**1.5*(x[b+i]**0.5\
                +x[b+i+1]**0.5)**3)/T0
            H[b+i+1,b+i+1]+= 2*xsi*(x[b+i]**0.5+3*x[b+i+1]**0.5)/(4*x[b+i+1]**1.5*(x[b+i+1]**0.5\
                +x[b+i]**0.5)**3)/T0
            H[b+i,b+i+1]+= 2*xsi/(2*(x[b+i]*x[b+i+1])**0.5*(x[b+i]**0.5+x[b+i+1]**0.5)**3)/T0
            H[b+i+1,b+i]+= 2*xsi/(2*(x[b+i]*x[b+i+1])**0.5*(x[b+i+1]**0.5+x[b+i]**0.5)**3)/T0
        return sparse.csr_matrix(H)/expansion_factor_ab**2
    return Hessian

















#creates bounds to b 
def create_b_bounds(n_discretization,n_wheels):
    
    length_wheels = 2*n_wheels*(n_discretization-1)
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([-np.inf]*(n_discretization-1))
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*length_wheels)
    ub.extend([np.inf]*(2*n_discretization-1+length_wheels))
    bounds = Bounds(lb,ub)
    return bounds



















#defines equality constraint matrix F
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal 
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), 
# number of discretization
#Output constraint Matrix F
def create_Linear_Constraints(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #flattened vector coordinates
    b = n_discretization-1
    u=n_discretization+n_discretization-1
    
    
    F = sparse.lil_matrix((4*(n_discretization-1),(2*n_wheels+1)*(n_discretization-1)+n_discretization))
    lb = np.zeros(4*(n_discretization-1))
    ub = np.zeros(4*(n_discretization-1))
    
    
    #iterate over each section to have the dynamics constraint on the 
    #two cartesian coordinates 
    #and the differential constraint
    for i in range(n_discretization-1):
        
        #first Cartesian coordinate dynamic constraint
        F[i,i]=-M_t[i][0]/expansion_factor_ab
        F[i,b+i]=-C_t[i][0]/2/expansion_factor_ab
        F[i,b+i+1]=-C_t[i][0]/2/expansion_factor_ab
        for j in range(2*n_wheels):
            F[i,u+i+j*(n_discretization-1)]=R_t[i][0][j]/expansion_factor_u  
        
        
        #second Cartesian coordinate dynamic constraint
        F[n_discretization-1+i,i]=-M_t[i][1]/expansion_factor_ab
        F[n_discretization-1+i,b+i]=-C_t[i][1]/2/expansion_factor_ab
        F[n_discretization-1+i,b+i+1]=-C_t[i][1]/2/expansion_factor_ab
        for j in range(2*n_wheels):
            F[n_discretization-1+i,u+i+j*(n_discretization-1)]=R_t[i][1][j]/expansion_factor_u  
            
            
            
        #Third Cartesian coordinate dynamic constraint
        F[2*(n_discretization-1)+i,i]=-M_t[i][2]/expansion_factor_ab
        F[2*(n_discretization-1)+i,b+i]=-C_t[i][2]/2/expansion_factor_ab
        F[2*(n_discretization-1)+i,b+i+1]=-C_t[i][2]/2/expansion_factor_ab
        for j in range(2*n_wheels):
            F[2*(n_discretization-1)+i,u+i+j*(n_discretization-1)]=R_t[i][2][j]/expansion_factor_u    

        
        
        #Differential contraint
        F[3*(n_discretization-1)+i,i]=2*1/(n_discretization-1)/expansion_factor_ab
        F[3*(n_discretization-1)+i,b+i]=1/expansion_factor_ab
        F[3*(n_discretization-1)+i,b+1+i]=-1/expansion_factor_ab
    return LinearConstraint(sparse.csr_matrix(F),lb,ub)











#defines innequality constraint (friction circle)
#Input friction coef mu, mass of the vehicle m, number of discretizations
#Output 1d vector remainder, which goes to zero when the inequality holds
def create_friction_circle(mu,mass,n_discretization,n_wheels,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization+n_discretization-1
    
    lb = np.full(n_wheels*(n_discretization-1),-np.inf)
    ub = np.zeros(n_wheels*(n_discretization-1))

    def friction_circle(t):
         #expand space to facilitate the solver
         
        decision_variables = t/expansion_factor_u
        remainder = np.zeros(n_wheels*(n_discretization-1))
        
        #in each section
        for i in range(n_discretization-1):
            
            #for every wheel
            for j in range(n_wheels):
                remainder[i+j*(n_discretization-1)] = (\
                    decision_variables[u+i+2*j*(n_discretization-1)]**2\
                        +decision_variables[u+i+(2*j+1)*(n_discretization-1)]
                             **2)-(mu*mass*9.81/n_wheels)**2
            
        return remainder
    return friction_circle, lb, ub
















#creates friction circle constraints gradient
def create_friction_circle_jacobian(n_discretization,n_wheels,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization+n_discretization-1
    
    def friction_circle_jac(t):
        #expand space to facilitate the solver
        x = t/expansion_factor_u
        B1=sparse.lil_matrix((n_wheels*(n_discretization-1),\
            (2*n_wheels+1)*(n_discretization-1)+n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            for j in range(n_wheels):
                B1[i+j*(n_discretization-1),u+i+2*j*(n_discretization-1)] = 2*x[u+i+2*j*(n_discretization-1)]
                B1[i+j*(n_discretization-1),u+i+(2*j+1)*(n_discretization-1)] = 2*x[u+i+(2*j+1)*(n_discretization-1)]
        return sparse.csr_matrix(B1)/expansion_factor_u
    return friction_circle_jac








#creates friction circle constraints gradient
def create_friction_circle_Hessian(n_discretization,n_wheels,expansion_factor_u):
    
    #flattened vector coordinates
    u=n_discretization+n_discretization-1
    
    def friction_circle_Hessian(x,v):
        
        B=sparse.lil_matrix(((2*n_wheels+1)*(n_discretization-1)+n_discretization,\
                (2*n_wheels+1)*(n_discretization-1)+n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            for j in range(n_wheels):
                Hessian_i = sparse.lil_matrix(((2*n_wheels+1)*(n_discretization-1)+n_discretization,\
                (2*n_wheels+1)*(n_discretization-1)+n_discretization))
                Hessian_i[u+i+2*j*(n_discretization-1),u+i+2*j*(n_discretization-1)]=2
                Hessian_i[u+i+(2*j+1)*(n_discretization-1),u+i+(2*j+1)*(n_discretization-1)]=2

                B+=v[i]*Hessian_i
        return sparse.csr_matrix(B)/expansion_factor_u**2
    return friction_circle_Hessian







#Helps building innitial guess of constant b and normalization T0, E0
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal 
#C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0,R_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u):
    
    #creates innitial guess
    x0 = np.zeros(n_discretization-1)
    x0 = np.append(x0,(np.ones(n_discretization))*b0*expansion_factor_ab)
    x0 = np.append(x0,np.zeros((2*n_wheels)*(n_discretization-1)))
    
    #flattened vector coordinates first 3 forces
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    u3=u2+n_discretization-1
    
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        u = np.linalg.pinv(R_t[i])@C_t[i]*b0*expansion_factor_u
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]
        x0[u3+i]=u[2]
        
    return x0




















#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output result
def optimization_SQP_abu_3(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):


    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    expansion_factor_ab = 1E3
    expansion_factor_u = 1E0
    
    
    #creating objective vector
    T0=1
    E0=1
    

    #create the objective information
    obj = create_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    obj_grad = create_gradient_objective(xsi,A_t,T0,E0,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    obj_hess = create_Hessian_objective(xsi,T0,n_discretization,n_wheels,expansion_factor_ab)


    #create bounds
    bounds = create_b_bounds(n_discretization,n_wheels)


    #Create Linear Constraints
    Linear_c = create_Linear_Constraints(R_t,M_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)
    
    #create Firction Circle Constraints
    friction_circle, lb_fc, ub_fc = create_friction_circle(mu,mass,n_discretization,n_wheels,expansion_factor_u)
    Grad_fc = create_friction_circle_jacobian(n_discretization,n_wheels,expansion_factor_u)
    Hessian_fc = create_friction_circle_Hessian(n_discretization,n_wheels,expansion_factor_u)
    Non_linear_c = NonlinearConstraint(friction_circle,lb_fc,ub_fc,\
                    Grad_fc, hess=Hessian_fc)



    
    b0=1e-4
    x0=build_x0(b0,R_t,C_t,n_discretization,n_wheels,expansion_factor_ab,expansion_factor_u)


    options = {
    'verbose': display,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
        }   
    

    result = scp.optimize.minimize(obj, x0, method='trust-constr',jac = obj_grad, hess = obj_hess,
         constraints=[Linear_c, Non_linear_c],options=options,bounds=bounds)
    
    decision_variables = result.x
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (friction_circle(decision_variables)<= 1E-6).all())
        print("Test friction circle initial guess ", (friction_circle(x0)<= 1E-6).all())
        
        
        #check the rescaling if you need the variables other than b
    return  decision_variables/expansion_factor_ab
