import numpy as np
import cvxpy as cp
from scipy import sparse 


#defines the objective vector
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio 
#vector A_t),number of discretization
#Output objective vector f
def create_objective_vector(xsi, A_t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    d = u2+n_discretization-1+n_discretization
    
    f = np.zeros(4*(n_discretization-1)+2*n_discretization)
    
    
    #for each path section
    for i in range(n_discretization-1):
        f[u1+i]=(1-xsi)*A_t[i][0]/E0
        f[u2+i]=(1-xsi)*A_t[i][1]/E0
        f[d+i]=2*xsi/T0
    return f













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
    
    
    F = sparse.lil_matrix((3*(n_discretization-1),2*n_discretization+4*
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
    return sparse.csr_matrix(F)













#creates bounds to b 
def create_b_bounds(x,n_discretization):
    #flattened vector coordinates
    b = n_discretization-1
    
    #create soc constraint vector
    soc_constraints = []
    #create all the b>=0 constraint
    for i in range(n_discretization):
        c_vec=np.zeros(2*n_discretization+4*(n_discretization-1))
        c_vec[b+i] = 1
        
        
        soc_constraints.append(cp.SOC(c_vec.T@x, cp.Constant(np.zeros(2))))
        
    return soc_constraints














#creates b/c cone constraint
def create_b_c_cones(x,n_discretization):
    
    #flattened vector coordinates
    b = n_discretization-1
    c = b+n_discretization+2*(n_discretization-1)
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the b/c constraint
    for i in range(n_discretization):
        
        #build the cone vector c_vec, which is 1 for b_k and 0 otherwise
        c_vec=np.zeros(2*n_discretization+4*(n_discretization-1))
        c_vec[b+i] = 1
        
        
        #build the cone matrix A_matrix, which is 2 for c_k in the first 
        # line, 1 for b_k in thensecond line, and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+4*(n_discretization-1)))
        A_matrix[0,c+i]=2
        A_matrix[1,b+i]=1
        
        
        #build the cone vector b_vec, which is -1 on the second line and 0 otherwise
        b_vec = np.zeros(2)
        b_vec[1] = -1
        
        
        soc_constraints.append(cp.SOC(c_vec.T @ x+cp.Constant(1), sparse.csr_matrix(A_matrix)@x+b_vec))
        
    return soc_constraints















#creates c/d cone constraint
def create_c_d_cones(x,n_discretization):
    
    #flattened vector coordinates
    c = n_discretization+3*(n_discretization-1)
    d=c+n_discretization
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the c/d constraint
    for i in range(n_discretization-1):
        
        #build the cone vector c_vec, which is 1 for d_k, c_{k+1}, c_k 
        #and 0 otherwise
        c_vec=np.zeros(2*n_discretization+4*(n_discretization-1))
        c_vec[d+i]=1
        c_vec[c+i]=1
        c_vec[c+1+i]=1
        
        
        #build the cone matrix A_matrix, which is 1 for d_k, -1 for c_{k+1}
        #, -1 for c_k in the
        #second line, and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+4*(n_discretization-1)))
        A_matrix[1,d+i]=1
        A_matrix[1,c+i]=-1
        A_matrix[1,c+1+i]=-1
        
        #build the cone vector b_vec, which is 2 on the first line and 
        #0 otherwise
        b_vec = np.zeros(2)
        b_vec[0] = 2
        
        
        soc_constraints.append(cp.SOC(c_vec.T @ x,sparse.csr_matrix(A_matrix)@x+cp.Constant(b_vec)))
        
    return soc_constraints

















#creates friction circle constraints
def create_friction_circle_cones(x,n_discretization,m,mu):
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the frisction circle constraints
    for i in range(n_discretization-1):
        
        #build the cone matrix A_matrix, which is 1 for u_1k on the first 
        #line and for u_2k on the second line,and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+4*(n_discretization-1)))
        A_matrix[0,u1+i]=1
        A_matrix[1,u2+i]=1
        
        
        soc_constraints.append(cp.SOC(cp.Constant(m*mu*9.81), sparse.csr_matrix(A_matrix)@x))
        
    return soc_constraints















#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_SOCP_abu(R_t,M_t,C_t,A_t,n_discretization,xsi,n_wheels,display):
    
    #create the decision variables vector
    x = cp.Variable(2*n_discretization+4*(n_discretization-1))
    
    
    #creating objective vector
    T0=1
    E0=1
    f = create_objective_vector(xsi, A_t,T0,E0,n_discretization)
    
    
    #creating equality constraint variables
    F=create_equality_constraint_matrix(R_t,M_t,C_t,n_discretization)
    g = np.zeros(3*(n_discretization-1))
    
    
    #creating cone constraints
    soc_constraints = []
  
    soc_constraints.extend(create_b_bounds(x,n_discretization))
    soc_constraints.extend(create_b_c_cones(x,n_discretization))
    soc_constraints.extend(create_c_d_cones(x,n_discretization))
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    soc_constraints.extend(create_friction_circle_cones(x,n_discretization,mass,mu))
    
    
    #set the SOCP problem
    prob = cp.Problem(cp.Minimize(f.T@x),soc_constraints+[F @ x == g])
    prob.solve(solver=cp.CLARABEL)

    # Print result.
    if display:
        print("Optimization terminated successfully")
        print(f"The optimal value is, {prob.value:.4f}")
    return  x.value