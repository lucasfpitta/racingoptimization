import numpy as np
import scipy as scp
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import BFGS
from scipy import sparse 














#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    b=n_discretization-1
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    def objective_function(decision_variables):
        
        cost=0
        
        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+(2*xsi/((decision_variables[b+i+1]**0.5+decision_variables[b+i]
                    **0.5)*T0)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+
                                        decision_variables[u2+i]*A_t[i][1])/E0)
        return cost
    
    return objective_function














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

        for i in range(n_discretization-1):
            f[b+i] += 2*xsi*(-1/(2*np.sqrt(x[b+i])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0
            f[b+i+1] += 2*xsi*(-1/(2*np.sqrt(x[b+i+1])*(np.sqrt(x[b+i+1])+np.sqrt(x[b+i]))**2))/T0
        return f
    return Grad






#defines Hessian of the function at a linearization point
def create_Hessian_objective(xsi,T0,n_discretization):
    
    #flattened vector coordinates
    b=n_discretization-1
    
    def Hessian(x):
    
        H = sparse.lil_matrix((3*(n_discretization-1)+n_discretization,3*(n_discretization-1)+n_discretization))
        
        
        #for each path section
        for i in range(n_discretization-1):
            H[b+i,b+i]+=2*xsi*(x[b+i+1]**0.5+3*x[b+i]**0.5)/(4*x[b+i]**1.5*(x[b+i]**0.5\
                +x[b+i+1]**0.5)**3)/T0
            H[b+i+1,b+i+1]+= 2*xsi*(x[b+i]**0.5+3*x[b+i+1]**0.5)/(4*x[b+i+1]**1.5*(x[b+i+1]**0.5\
                +x[b+i]**0.5)**3)/T0
            H[b+i,b+i+1]+= 2*xsi/(2*(x[b+i]*x[b+i+1])**0.5*(x[b+i]**0.5+x[b+i+1]**0.5)**3)/T0
            H[b+i+1,b+i]+= 2*xsi/(2*(x[b+i]*x[b+i+1])**0.5*(x[b+i+1]**0.5+x[b+i]**0.5)**3)/T0
        return sparse.csr_matrix(H)
    return Hessian

















#creates bounds to b 
def create_b_bounds(n_discretization):
    
    lb=[]
    ub=[]
    
    #lower bounds above 0 to avoid objective problems
    lb.extend([-np.inf]*(n_discretization-1))
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*2*(n_discretization-1))
    ub.extend([np.inf]*(4*n_discretization-3))
    bounds = Bounds(lb,ub)
    return bounds

















#defines equality constraint matrix F
#Input Force R_t (3d array with n_discretizatio matrix R_t), Mass and Centrifugal 
#M_t, C_t (2d array with n_discretizatio of vectors M_t and C_t), 
# number of discretization
#Output constraint Matrix F
def create_Linear_Constraints(R_t,M_t,C_t,n_discretization):
    
    #flattened vector coordinates
    b = n_discretization-1
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    F = sparse.lil_matrix((3*(n_discretization-1),n_discretization+3*
                  (n_discretization-1)))
    lb = np.zeros(3*(n_discretization-1))
    ub = np.zeros(3*(n_discretization-1))
    
    
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
    return LinearConstraint(sparse.csr_matrix(F),lb,ub)











def create_friction_circle(mu,mass,n_discretization): 
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    lb = np.full(n_discretization-1,-np.inf)
    ub = np.zeros(n_discretization-1)
    
    def friction_circle(x):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = (x[u1+i]**2+x[u2+i]
                             **2)**0.5-mu*mass*9.81
        return remainder
    return friction_circle, lb, ub
















#creates friction circle constraints gradient
def create_friction_circle_jacobian(n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=u1+n_discretization-1
    
    def friction_circle_jac(x):
        B1=sparse.lil_matrix((n_discretization-1,n_discretization+3*(n_discretization-1)))
        
        #create all the frisction circle constraints
        for i in range(n_discretization-1):
            B1[i,u1+i] = x[u1+i]/(x[u1+i]**2+x[u2+i]**2)**0.5
            B1[i,u2+i] = x[u2+i]/(x[u1+i]**2+x[u2+i]**2)**0.5
        return sparse.csr_matrix(B1)
    return friction_circle_jac
















#Helps building innitial guess of constant b and normalization T0, E0
#Input Force R_t (3d array with n_discretizatio matrix R_t), Centrifugal 
#C_t (2d array with n_discretizatio of vector C_t), number of discretization
#Output 1d flattened vector of initial guess
def build_x0(b0,R_t,C_t,n_discretization):
    
    #creates innitial guess
    x0 = np.zeros(n_discretization-1)
    x0 = np.append(x0,np.ones(n_discretization)*b0)
    x0 = np.append(x0,np.zeros(2*(n_discretization-1)))

    
    #flattened vector coordinates
    b=n_discretization-1
    u1=b+n_discretization
    u2=u1+n_discretization-1
    
    
    #calculates forces that are necessary for constant u
    for i in range(n_discretization-1):
        u = (x0[b+i+1]+x0[b+i])/2*\
            np.linalg.inv(R_t[i])@C_t[i]
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]

    return x0





















#Optimizer
#Input Force R_t (3d array with n_discretizatio matrix R_t), Power, Mass and 
#Centrifugal A_t, M_t, C_t (2d array with n_discretizatio of vectors A_t, 
#M_t and C_t), number of discretization, xsi optimization scalar
#Output result
def optimization_SQP_abu(R_t,M_t,C_t,A_t,n_discretization,xsi,display):


    mu=1 #friction coeficient
    mass=85 #mass of the vehicle
    epsilon=1e-3 #hypercube size
    
    
    #creating objective vector
    T0=1
    E0=1
    

    #create the objective information
    obj = create_objective(xsi,A_t,T0,E0,n_discretization)
    obj_grad = create_gradient_objective(xsi,A_t,T0,E0,n_discretization)
    obj_hess = create_Hessian_objective(xsi,T0,n_discretization)


    #create bounds
    bounds = create_b_bounds(n_discretization)


    #Create Linear Constraints
    Linear_c = create_Linear_Constraints(R_t,M_t,C_t,n_discretization)
    
    #create Firction Circle Constraints
    friction_circle, lb_fc, ub_fc = create_friction_circle(mu,mass,n_discretization)
    Grad_fc = create_friction_circle_jacobian(n_discretization)
    #Hessian_fc = create_Hessian_objective(xsi,T0,n_discretization)
    Non_linear_c = NonlinearConstraint(friction_circle,lb_fc,ub_fc,\
                    Grad_fc, hess=BFGS())



    
    b0=1e-4
    x0=build_x0(b0,R_t,C_t,n_discretization)


    options = {
    'verbose': True,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
        }   
    

    result = scp.optimize.minimize(obj, x0, method='trust-constr',jac = obj_grad, hess = obj_hess,
         constraints=[Linear_c, Non_linear_c],options=options,bounds=bounds)
    
    decision_variables = result.x
    if display:
        print("T0 ", T0, " E0 ", E0)
        print("Test friction circle ", (friction_circle(decision_variables)<= 1E-6).all())
        print("Test friction circle initial guess ", (friction_circle(x0)<= 1E-6).all())
    return  decision_variables














































"""

def create_constraint3(mu,mass,n_discretization):
    
    #flattened vector coordinates
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    
    def constraint3(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = (decision_variables[u1+i]**2+decision_variables[u2+i]
                             **2)**0.5-mu*mass*9.81
            
        return remainder
    return constraint3











def create_Hessian_constraint(n_discretization):
    
    #flattened vector coordinates
    u1 = 2*n_discretization-1
    u2 = u1+n_discretization-1
    mu = n_discretization+6*(n_discretization-1)
    
    def Hessian_constraint(x):
    
        H = np.zeros((3*(n_discretization-1)+n_discretization,3*(n_discretization-1)+n_discretization))
        
        
        #for each path section
        for i in range(n_discretization-1):
            H[u1+i][u1+i]=x[mu+i]*x[u2+i]**2/(x[u1+i]**2+x[u2+i]**2)**1.5
            H[u1+i][u2+i]=-x[mu+i]*x[u1+i]*x[u2+i]/(x[u1+i]**2+x[u2+i]**2)**1.5
            H[u2+i][u2+i]=x[mu+i]*x[u1+i]**2/(x[u1+i]**2+x[u2+i]**2)**1.5
            H[u2+i][u1+i]=-x[mu+i]*x[u1+i]*x[u2+i]/(x[u1+i]**2+x[u2+i]**2)**1.5
        return H
    return Hessian_constraint









#defines the objective
#Input optimization weight scalar xsi, Power A_t (2d array with 
#n_discretizatio vector A_t), 
# number of discretization
#Output cost of the solution (scalar)
def create_objective(xsi,A_t,T0,E0,n_discretization):
    
    #flattened vector coordinates
    b=n_discretization-1
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    
    def objective_function(decision_variables):
        
        cost=0
        
        #sum over the path 
        for i in range(n_discretization-1):
            cost = cost+(2*xsi/((decision_variables[b+i+1]**0.5+decision_variables[b+i]
                    **0.5)*T0)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+
                                        decision_variables[u2+i]*A_t[i][1])/E0)
        return cost
    
    return objective_function























def teste_Grad(grad_anal,f,epsilon,x):
    grad_num = np.zeros(len(grad_anal))
    for i in range (len(grad_anal)):
        x_up = x.copy()
        x_up[i]=x_up[i]+epsilon
        x_down = x.copy()
        x_down[i]=x_down[i]-epsilon
        f_up = f(x_up)
        f_down = f(x_down)
        grad_num[i]=1/(2*epsilon)*(f_up-f_down)

    return grad_num






def teste_Hessian(hess_anal,f,epsilon,x):
    hess_num = np.zeros(np.shape(hess_anal))
    for i in range (np.shape(hess_anal)[0]):
        for j in range (np.shape(hess_anal)[1]):
            e_i = np.zeros(np.shape(hess_anal)[0])
            e_j = np.zeros(np.shape(hess_anal)[1])
            e_i[i] = epsilon
            e_j[j] = epsilon

            f_pp = f(x+e_i+e_j)
            f_pm = f(x+e_i-e_j)
            f_mp = f(x-e_i+e_j)
            f_mm = f(x-e_i-e_j)

            hess_num[i][j]=1/(4*epsilon**2)*(f_pp-f_pm-f_mp+f_mm)
    
    output_file = "matrix_export.csv"
    np.savetxt(output_file, hess_anal-hess_num, delimiter=",", fmt="%g")
    print()
    print("Numerico",np.max(np.abs(hess_num)))
    print()
    print("Analitico",np.max(np.abs(hess_anal)))
    print()
    print("Erro",np.max(np.abs(hess_anal-hess_num)))
    print()

    return np.max(np.abs(hess_anal-hess_num)/np.max(np.abs(hess_anal)))










def teste_Grad_friction(grad_anal,f,epsilon,x):
    grad_num = np.zeros(np.shape(grad_anal))
    for j in range(np.shape(grad_anal)[0]):
        for i in range (np.shape(grad_anal)[1]):
            x_up = x.copy()
            x_up[i]=x_up[i]+epsilon
            x_down = x.copy()
            x_down[i]=x_down[i]-epsilon
            f_up = f(x_up)[j]
            f_down = f(x_down)[j]
            grad_num[j][i]=1/(2*epsilon)*(f_up-f_down)

    print()
    print("Numerico",np.max(np.abs(grad_num)))
    print()
    print("Analitico",np.max(np.abs(grad_anal)))
    print()
    print("Erro",np.max(np.abs(grad_anal-grad_num)))
    print()

    return np.max(np.abs(grad_anal-grad_num))

    
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
    Hessian_c = create_Hessian_constraint(n_discretization)
    B1 = create_friction_circle(n_discretization)

    


    # obj = create_objective(xsi,A_t,T0,E0,n_discretization)
    epsilon =1e-5
    # for i in range(10):
    #     x_test = np.random.rand(n_discretization+3*(n_discretization-1))
    #     grad_anal = Grad_f(x_test)
    #     print("teste Grad",teste_Grad(grad_anal,obj,epsilon,x_test)<=1e-5)

    
    # for i in range(100):
    #     x_test = np.random.rand(n_discretization+3*(n_discretization-1))
    #     hess_anal = Hessian(x_test)
    #     erros [i] = teste_Hessian(hess_anal,obj,epsilon,x_test)
        
        
    erros = np.zeros(100)
    fc = create_constraint3(mu,mass,n_discretization)
    for i in range(100):
        x_test = np.random.rand(n_discretization+3*(n_discretization-1))
        grad_anal = B1(x_test)
        print("teste Fric",teste_Grad_friction(grad_anal,fc,epsilon,x_test)<=1e-5)
        erros [i] = teste_Grad_friction(grad_anal,fc,epsilon,x_test)
        
    print("Erro relativo máximo: ", np.max(erros))
    print("Erro relativo médio: ", np.average(erros))


    
    #creates the 0 blocks
    Block1 = np.zeros((3*(n_discretization-1),3*(n_discretization-1)))
    Block2 = np.zeros((3*(n_discretization-1),2*n_discretization-1)) 
    Block3 = np.zeros((2*n_discretization-1,3*(n_discretization-1)))
    Block4 = np.zeros((2*n_discretization-1,2*n_discretization-1)) 
    
    d=0
    while np.linalg.norm(deltaX1-deltaX0)>1E-9 or d>=1000:
        print(f"iteration {d}")

        for i in range(n_discretization):
            if deltaX1[n_discretization-1+i]<=0:
                deltaX1[n_discretization-1+i]=1e-10

        Hess = Hessian(deltaX1[0:n_discretization+3*(n_discretization-1)])
        Hess_c = Hessian_c(deltaX1)
        B_h= np.vstack((B0,B1(deltaX1[0:n_discretization+3*(n_discretization-1)])))
        B_g= np.vstack((B0,-mu*mass*9.81+B1(deltaX1[0:n_discretization+3*(n_discretization-1)])))
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
            [Hess+Hess_c, np.transpose(F), np.transpose(B_h)],
            [F, Block1, Block2],
            [B_h, Block3, Block4]
            ])
        



        deltaX0 = deltaX1.copy()
        deltaX1 = deltaX0-(scp.linalg.lu(Lag_Hessian)[0]@Lag_Grad).reshape(-1)
        d=d+1
        


    # Print result.
    if display:
        print("Optimization terminated successfully")
    return  np.ones(n_discretization)



















"""