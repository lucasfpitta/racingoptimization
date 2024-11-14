import numpy as np
import scipy as scp
import cvxpy as cp



#defines gradient of the function at a linearization point
def create_gradient_objective(xsi,A_t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    b=n_discretization-1
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    
    def Grad(x):
        f = np.zeros(3*(n_discretization-1)+n_discretization)
        
        #for each path section
        for i in range(n_discretization-1):
            f[u1+i]=(1-xsi)*A_t[i][0]/E0
            f[u2+i]=(1-xsi)*A_t[i][1]/E0
            f[b+i]=f[b+i]+2*xsi/T0*(-1/(2*x[b+i]**0.5*(x[b+i]**0.5+x[b+1+i]**0.5)**2))
            f[b+i+1]=2*xsi/T0*(-1/(2*x[b+i+1]**0.5*(x[b+i+1]**0.5+x[b+i]**0.5)**2))
        return f
    return Grad






#defines Hessian of the function at a linearization point
def create_Hessian_objective(xsi,T0,n_discretization):
    
    #flattened vector coordinates
    b=n_discretization-1
    
    def Hessian(x):
    
        H = np.zeros((3*(n_discretization-1)+n_discretization,3*(n_discretization-1)+n_discretization))
        
        
        #for each path section
        for i in range(n_discretization-1):
            H[b+i][b+i]=H[b+i][b+i]+2*xsi/T0*(x[b+i+1]**0.5+3*x[b+i]**0.5)/(4*x[b+i]**1.5*(x[b+i]**0.5\
                +x[b+i+1]**0.5)**3)
            H[b+i+1][b+i+1]= 2*xsi/T0*(x[b+i]**0.5+3*x[b+i+1]**0.5)/(4*x[b+i+1]**1.5*(x[b+i+1]**0.5\
                +x[b+i]**0.5)**3)
            H[b+i][b+i+1] = 2*xsi/T0/(2*(x[b+i]*x[b+i+1])**0.5*(x[b+i]**0.5+x[b+i+1]**0.5)**3)
            H[b+i+1][b+i] = 2*xsi/T0/(2*(x[b+i]*x[b+i])**0.5*(x[b+i+1]**0.5+x[b+i]**0.5)**3)
        return H
    return Hessian





#defines equality constraint matrix F
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal 
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), 
# number of discretization
#Output constraint Matrix F
def create_equality_constraint_matrix(R_t,M_t,C_t,n_discretization):
    
    #flattened vector coordinates
    b = n_discretization-1
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    F = np.zeros((3*(n_discretization-1),n_discretization+3*
                  (n_discretization-1)))
    
    
    #iterate over each section to have the dynamics constraint on the 
    #two cartesian coordinates 
    #and the differential constraint
    for i in range(n_discretization-1):
        
        #first Cartesian coordinate dynamic constraint
        F[i,i]=-M_t[i][0]
        F[i,b+i]=-C_t[i][0]/2
        F[i,b+i+1]=-C_t[i][0]/2
        F[i,u1+i]=R_t[i][0][0]
        F[i,u2+i]=R_t[i][0][1]
        
        
        #second Cartesian coordinate dynamic constraint
        F[n_discretization-1+i,i]=-M_t[i][1]
        F[n_discretization-1+i,b+i]=-C_t[i][1]/2
        F[n_discretization-1+i,b+i+1]=-C_t[i][1]/2
        F[n_discretization-1+i,u1+i]=R_t[i][1][0]
        F[n_discretization-1+i,u2+i]=R_t[i][1][1]
        
        
        #Differential contraint
        F[2*(n_discretization-1)+i,i]=2*1/(n_discretization-1)
        F[2*(n_discretization-1)+i,b+i]=1
        F[2*(n_discretization-1)+i,b+1+i]=-1
    return F













#creates bounds to b 
def create_b_bounds(n_discretization):
    #flattened vector coordinates
    b = n_discretization-1
    
    #create soc constraint vector
    B0=np.zeros((n_discretization,n_discretization+3*(n_discretization-1)))
    #create all the b>=0 constraint
    for i in range(n_discretization):
        B0[i][b+i]=1
    return B0



















#creates friction circle constraints
def create_friction_circle(n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    
    def friction_circle(x):
        B1=np.zeros((n_discretization-1,n_discretization+3*(n_discretization-1)))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            B1[i][u1+i] = x[u1+i]/(x[u1+i]**2+x[u2+i]**2)**0.5
            B1[i][u2+i] = x[u2+i]/(x[u1+i]**2+x[u2+i]**2)**0.5
        return B1
    return friction_circle















#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output result
def optimization_SQP_abu(R_t,M_t,C_t,A_t,n_discretization,xsi,display):


    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    #create initial guess
    x0 = 1e-8*np.ones(n_discretization+3*(n_discretization-1))
    deltaX0 = np.zeros(n_discretization+3*(n_discretization-1)+3*(n_discretization-1)+2*n_discretization-1)
    deltaX1 = np.append(x0,np.zeros(3*(n_discretization-1)+2*n_discretization-1))
    
    #creating objective vector
    T0=1
    E0=1
    
    #create constant Matrices
    B0=create_b_bounds(n_discretization)
    F=create_equality_constraint_matrix(R_t,M_t,C_t,n_discretization)
    
    #create the approx Matrices
    Grad_f = create_gradient_objective(xsi,A_t,T0,E0,n_discretization)
    Hessian = create_Hessian_objective(xsi,T0,n_discretization)
    B1 = create_friction_circle(n_discretization)
    
    
    #creates the 0 blocks
    Block1 = np.zeros((3*(n_discretization-1),3*(n_discretization-1)))
    Block2 = np.zeros((3*(n_discretization-1),2*n_discretization-1)) 
    Block3 = np.zeros((2*n_discretization-1,3*(n_discretization-1)))
    Block4 = np.zeros((2*n_discretization-1,2*n_discretization-1)) 
    
    d=0
    while np.linalg.norm(deltaX1-deltaX0)>1E-9 or d>=1000:
        print(f"iteration {d}")
        print(deltaX1)

        for i in range(n_discretization):
            if deltaX1[n_discretization-1+i]<=0:
                deltaX1[n_discretization-1+i]=1e-10

        Hess = Hessian(deltaX1[0:n_discretization+3*(n_discretization-1)])
        B_h= np.vstack((B0,B1(deltaX1[0:n_discretization+3*(n_discretization-1)])))
        B_g= np.vstack((B0,mu*mass*9.81-B1(deltaX1[0:n_discretization+3*(n_discretization-1)])))
        Grad = Grad_f(deltaX1[0:n_discretization+3*(n_discretization-1)])


        #Build the Lagrangian Gradient
        primal_variables = \
        np.transpose(Grad)+np.transpose(F)@deltaX1[n_discretization+3*(n_discretization-1):\
                                                   n_discretization+6*(n_discretization-1)]+\
        np.transpose(B_h)@deltaX1[n_discretization+6*(n_discretization-1):\
                                2*n_discretization+7*(n_discretization-1)]#+Hess@deltaX1[0:n_discretization+3*(n_discretization-1)]
    

        
        
        dual_variables1 = F@deltaX1[0:n_discretization+3*(n_discretization-1)]

        dual_variables2=B_g@deltaX1[0:n_discretization+3*(n_discretization-1)]

        Lag_Grad = np.vstack((primal_variables[:, None] ,dual_variables1[:, None] ,dual_variables2[:, None] ))


        #Build the Lagrangian Hessian
        Lag_Hessian =  np.block([
            [Hess, np.transpose(F), np.transpose(B_h)],
            [F, Block1, Block2],
            [B_h, Block3, Block4]
            ])
        



        deltaX0 = deltaX1
        deltaX1 = deltaX0-(scp.linalg.lu(Lag_Hessian)[0]@Lag_Grad).reshape(-1)
        d=d+1
        


    # Print result.
    if display:
        print("Optimization terminated successfully")
    return  np.ones(n_discretization)