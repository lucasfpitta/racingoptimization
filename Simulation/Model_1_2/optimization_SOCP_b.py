import numpy as np
import cvxpy as cp
from scipy import sparse 




# Create the vectors F_it for both the objective and constraints
# Input R_t (2d matrix with n_discretization vector), M_t (2d array with 
#n_discretization vector M_t)
# C_tM_t (2d array with n_discretization vector C_t), n_discretization
# Output F_1t and F_2t (2d array with n_discretization vector F_1t and vector 
# F_2t)
def create_F_it(R_t,M_t,C_t,n_discretization):
    
    #create the n_dicretization vectors F_1t and F_2t
    F_1t, F_2t = np.zeros(np.shape(M_t)), np.zeros(np.shape(M_t))
    
    #loop to build the vectors
    for i in range(n_discretization-1):
        F_1t[i] = np.linalg.inv(R_t[i])@(M_t[i]/(2*1/(n_discretization-1))+C_t[i]/2)
        F_2t[i] = np.linalg.inv(R_t[i])@(-M_t[i]/(2*1/(n_discretization-1))+C_t[i]/2)
    return F_1t, F_2t


















#defines the objective vector
#Input optimization weight scalar xsi, Power A_t (2d array with n_discretizatio
# vector A_t), force vectors F1_t, F_2t, number of discretization
#Output objective vector f
def create_objective_vector(xsi, A_t, F_1t, F_2t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    d = 2*n_discretization
    f = np.zeros(2*n_discretization+n_discretization-1)
     
    #for each path section
    f[0]=(1-xsi)*(F_2t[0].T@A_t[0])/E0
    
    
    for i in range(n_discretization-2):
        f[i+1]=(1-xsi)*(F_1t[i].T@A_t[i]+F_2t[i+1].T@A_t[i+1])/E0
        f[d+i]=2*xsi/T0
        
    f[n_discretization-1]=(1-xsi)*(F_1t[n_discretization-2].T@
                                   A_t[n_discretization-2])/E0
    f[-1]=f[d]=2*xsi/T0
    return f














#creates bounds to b 
def create_b_bounds(x,n_discretization):   
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the b>=0 constraint
    for i in range(n_discretization):
        c_vec=np.zeros(2*n_discretization+n_discretization-1)
        c_vec[i] = 1
        soc_constraints.append(cp.SOC(c_vec.T@x, cp.Constant(np.zeros(2))))
    return soc_constraints















#creates b/c cone constraint
def create_b_c_cones(x,n_discretization):
    #flattened vector coordinates
    c = n_discretization
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the b/c constraint
    for i in range(n_discretization):
        
        #build the cone vector c_vec, which is 1 for b_k and 0 otherwise
        c_vec=np.zeros(2*n_discretization+n_discretization-1)
        c_vec[i] = 1
        
        #build the cone matrix A_matrix, which is 2 for c_k in the first line, 
        # 1 for b_k in the second line, and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+n_discretization-1))
        A_matrix[0,c+i]=2
        A_matrix[1,i]=1
        
        
        #build the cone vector b_vec, which is -1 on the second line and 0 otherwise
        b_vec = np.zeros(2)
        b_vec[1] = -1
        
        
        soc_constraints.append(cp.SOC(c_vec.T @ x+cp.Constant(1), sparse.csr_matrix(A_matrix)@x+b_vec))
        
    return soc_constraints














#creates c/d cone constraint
def create_c_d_cones(x,n_discretization):
    
    #flattened vector coordinates
    c = n_discretization
    d=c+n_discretization
    
    #create soc constraint vector
    soc_constraints = []
    
    #create all the c/d constraint
    for i in range(n_discretization-1):
        #build the cone vector c_vec, which is 1 for d_k, c_{k+1}, c_k and 
        #0 otherwise
        c_vec=np.zeros(2*n_discretization+n_discretization-1)
        c_vec[d+i]=1
        c_vec[c+i]=1
        c_vec[c+1+i]=1
        
        
        #build the cone matrix A_matrix, which is 1 for d_k, -1 for c_{k+1}, 
        # -1 for c_k in the second line, and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+n_discretization-1))
        A_matrix[1,d+i]=1
        A_matrix[1,c+i]=-1
        A_matrix[1,c+1+i]=-1
        
        #build the cone vector b_vec, which is 2 on the first line and 0 otherwise
        b_vec = np.zeros(2)
        b_vec[0] = 2
        
        soc_constraints.append(cp.SOC(c_vec.T @ x, sparse.csr_matrix(A_matrix)@x+cp.Constant(b_vec)))
        
    return soc_constraints













#creates friction circle constraints
def create_friction_circle_cones(x,F_1t,F_2t,n_discretization,m,mu):
        
    #create soc constraint vector
    soc_constraints = []
    
    #create all the frisction circle constraints
    for i in range(n_discretization-1):
        
        #build the cone matrix A_matrix, which is F2[k][0] for b_k and  
        # F1[k][0] for b_{k+1} on the first line and F2[k][1] for b_k and  
        # F1[k][1] for b_{k+1} second line, and 0 otherwise
        A_matrix = sparse.lil_matrix((2,2*n_discretization+n_discretization-1))
        A_matrix[0,i]=F_2t[i][0]
        A_matrix[0,i+1]=F_1t[i][0]
        A_matrix[1,i]=F_2t[i][1]
        A_matrix[1,i+1]=F_1t[i][1]
            
        soc_constraints.append(cp.SOC(cp.Constant(m*mu*9.81), sparse.csr_matrix(A_matrix)@x))
        
    return soc_constraints













#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output scipy result and innitial guess x0
def optimization_SOCP_b(R_t,M_t,C_t,A_t,n_discretization,xsi,display):
    
    #create the decision variables vector
    x = cp.Variable(2*n_discretization+n_discretization-1)
    
    #creating objective vector
    T0=1
    E0=1
    
    #Creating force matrices F_1t and F_2t
    F_1t, F_2t = create_F_it(R_t,M_t,C_t,n_discretization)
    
    f = create_objective_vector(xsi, A_t, F_1t, F_2t,T0,E0,n_discretization)
    
    #creating cone constraints
    soc_constraints = []
    soc_constraints.extend(create_b_bounds(x,n_discretization))
    soc_constraints.extend(create_b_c_cones(x,n_discretization))
    soc_constraints.extend(create_c_d_cones(x,n_discretization))
    
    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    
    soc_constraints.extend(create_friction_circle_cones(x,F_1t, F_2t,
                                            n_discretization,mass,mu))
    
    #set the SOCP problem
    prob = cp.Problem(cp.Minimize(f.T@x),soc_constraints)
    prob.solve(solver = cp.CLARABEL)

    # Print result.
    if display:
        print("Optimization terminated successfully")
        print(f"The optimal value is, {prob.value:.4f}")
    return  x.value