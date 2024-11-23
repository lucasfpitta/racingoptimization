import numpy as np
import scipy as scp
from scipy.optimize import Bounds
from scipy.optimize import NonlinearConstraint
from scipy import sparse 



# Create the vectors F_it for both the objective and constraints
# Input R_t (2d matrix with n_discretization vector), M_t (2d array 
# with n_discretization vector M_t)
# C_tM_t (2d array with n_discretization vector C_t), n_discretization
# Output F_1t and F_2t (2d array with n_discretization vector F_1t and vector F_2t)

def create_F_it(R_t,M_t,C_t,n_discretization):
    
    #create the n_dicretization vectors F_1t and F_2t
    F_1t, F_2t = np.zeros(np.shape(M_t)), np.zeros(np.shape(M_t))
    
    #loop to build the vectors
    for i in range(n_discretization-1):
        F_1t[i] = np.linalg.inv(R_t[i])@(M_t[i]/(2*1/(n_discretization-1))+
                                         C_t[i]/2)
        F_2t[i] = np.linalg.inv(R_t[i])@(-M_t[i]/(2*1/(n_discretization-1))+
                                         C_t[i]/2)
    return F_1t, F_2t














#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi, A_t, F_1t,F_2t,T0,E0,n_discretization,expansion_factor):    
    def objective_function(b):   
        #expand space to facilitate the solver
        decision_variables = b/expansion_factor
        cost=0

        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+2*xsi/(decision_variables[i+1]**0.5+decision_variables[i]**
                0.5)/T0+(1-xsi)*np.transpose(decision_variables[i+1]*F_1t[i]+
                                          decision_variables[i]*F_2t[i])@A_t[i]/E0
        return cost
    
    return objective_function














#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,F_1t,F_2t,T0,E0,n_discretization,expansion_factor):
    def Grad(b):
        
        x = b/expansion_factor
        f = np.zeros(n_discretization)

        for i in range(n_discretization-1):
            f[i] += 2*xsi*(-1/(2*np.sqrt(x[i])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0+\
                (1-xsi)/E0*np.transpose(F_2t[i])@A_t[i]
            f[i+1] += 2*xsi*(-1/(2*np.sqrt(x[i+1])*(np.sqrt(x[i+1])+np.sqrt(x[i]))**2))/T0+\
                (1-xsi)/E0*np.transpose(F_1t[i])@A_t[i]
        return f/expansion_factor
    return Grad






#defines Hessian of the function at a linearization point
def create_Hessian_objective(xsi,T0,n_discretization,expansion_factor):    
    def Hessian(b):
        x = b/expansion_factor
        H = sparse.lil_matrix((n_discretization,n_discretization))
        
        
        #for each path section
        for i in range(n_discretization-1):
            H[i,i]+=2*xsi*(x[i+1]**0.5+3*x[i]**0.5)/(4*x[i]**1.5*(x[i]**0.5\
                +x[i+1]**0.5)**3)/T0
            H[i+1,i+1]+= 2*xsi*(x[i]**0.5+3*x[i+1]**0.5)/(4*x[i+1]**1.5*(x[i+1]**0.5\
                +x[i]**0.5)**3)/T0
            H[i,i+1]+= 2*xsi/(2*(x[i]*x[i+1])**0.5*(x[i]**0.5+x[i+1]**0.5)**3)/T0
            H[i+1,i]+= 2*xsi/(2*(x[i]*x[i+1])**0.5*(x[i+1]**0.5+x[i]**0.5)**3)/T0
        return sparse.csr_matrix(H)/expansion_factor**2
    return Hessian

















#creates bounds to b 
def create_b_bounds(n_discretization):
    
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([1E-6]*n_discretization)
    ub.extend([np.inf]*n_discretization)
    bounds = Bounds(lb,ub)
    return bounds









def create_friction_circle(F_1t,F_2t,mu,mass,n_discretization,expansion_factor): 
    lb = np.full(n_discretization-1,-np.inf)
    ub = np.zeros(n_discretization-1)
    
    def friction_circle(b):
        x = b/expansion_factor
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = -(mu*mass*9.81)**2+np.dot((x[i+1]\
                    *F_1t[i]+x[i]*F_2t[i]),(x[i+1]\
                    *F_1t[i]+x[i]*F_2t[i]))
        return remainder
    return friction_circle, lb, ub
















#creates friction circle constraints gradient
def create_friction_circle_jacobian(F_1t,F_2t,n_discretization,expansion_factor):
    
    def friction_circle_jac(b):
        x = b/expansion_factor
        B1=sparse.lil_matrix((n_discretization-1,n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            B1[i,i] = 2*np.dot(F_2t[i],(x[i+1]
                    *F_1t[i]+x[i]*F_2t[i]))
            B1[i,i+1] = 2*np.dot(F_1t[i],(x[i+1]
                    *F_1t[i]+x[i]*F_2t[i])) 
        return sparse.csr_matrix(B1)/expansion_factor
    return friction_circle_jac








#creates friction circle constraints gradient
def create_friction_circle_Hessian(F_1t,F_2t,n_discretization,expansion_factor):
    def friction_circle_Hessian(x,v):
        B=sparse.lil_matrix((n_discretization,n_discretization))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            Hessian_i = sparse.lil_matrix((n_discretization,n_discretization))
            Hessian_i[i,i]=2*np.dot(F_2t[i],F_2t[i])
            Hessian_i[i+1,i]=2*np.dot(F_1t[i],F_2t[i])
            Hessian_i[i+1,i+1]=2*np.dot(F_1t[i],F_1t[i])
            Hessian_i[i,i+1]=2*np.dot(F_2t[i],F_1t[i])
            B+=v[i]*Hessian_i
        return sparse.csr_matrix(B)/expansion_factor**2
    return friction_circle_Hessian
























#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output result
def optimization_SQP_b(R_t,M_t,C_t,A_t,n_discretization,xsi,display):


    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    expansion_factor = 1E3
    
    
    #Creating force matrices F_1t and F_2t
    F_1t, F_2t = create_F_it(R_t,M_t,C_t,n_discretization)

    
    
    #creating objective vector
    T0=1
    E0=1
    

    #create the objective information
    obj = create_objective(xsi, A_t, F_1t,F_2t,T0,E0,n_discretization,expansion_factor)
    obj_grad = create_gradient_objective(xsi,A_t,F_1t,F_2t,T0,E0,n_discretization,expansion_factor)
    obj_hess = create_Hessian_objective(xsi,T0,n_discretization,expansion_factor)


    #create bounds
    bounds = create_b_bounds(n_discretization)

    
    #create Firction Circle Constraints
    friction_circle, lb_fc, ub_fc = create_friction_circle(F_1t,F_2t,mu,mass,n_discretization,expansion_factor)
    Grad_fc = create_friction_circle_jacobian(F_1t,F_2t,n_discretization,expansion_factor)
    Hessian_fc = create_friction_circle_Hessian(F_1t,F_2t,n_discretization,expansion_factor)
    Non_linear_c = NonlinearConstraint(friction_circle,lb_fc,ub_fc,\
                    Grad_fc, hess=Hessian_fc)



    
    b0=1e-4
    #creates initial guess inside the friction circle 
    x0 = np.ones(n_discretization)*expansion_factor*b0


    options = {
    'verbose': display,      # Display iteration info
    'xtol': 1e-6, 'gtol': 1e-6, 'maxiter': 1000,
    "initial_tr_radius": 1.0
        }   
    

    

    result = scp.optimize.minimize(obj, x0, method='trust-constr',jac = obj_grad, hess = obj_hess,
         constraints=Non_linear_c,options=options,bounds=bounds)
    
    decision_variables = result.x
    
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (friction_circle(decision_variables)<= 1E-6).all())
        print("Test friction circle initial guess ", (friction_circle(x0)<= 1E-6).all())

    return  decision_variables/expansion_factor
