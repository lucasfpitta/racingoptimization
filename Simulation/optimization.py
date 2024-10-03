import numpy as np
import scipy as scp

def create_objective(xsi, A_t,n_discretization):
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def objective_function(decision_variables):
        cost=0
        for i in range(n_discretization-1):
            cost = cost+2*xsi/(decision_variables[i+1]**0.5+decision_variables[i]**0.5)+(1-xsi)*(decision_variables[u1+i]*A_t[i][0]+decision_variables[u2+i]*A_t[i][1])
        return cost
    return objective_function

def create_constraint1(R_t,M_t,C_t,n_discretization):
    a = n_discretization
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def constraint1(decision_variables):
        remainder = np.zeros(2*(n_discretization-1))
        for i in range(n_discretization-1):
            remainder[i]=(R_t[i][0][0]*decision_variables[u1+i]+R_t[i][0][1]*decision_variables[u2+i])-M_t[i][0]*decision_variables[a+i]-C_t[i][0]*(decision_variables[i]+decision_variables[i+1])/2
            remainder[n_discretization-1+i]=(R_t[i][1][0]*decision_variables[u1+i]+R_t[i][1][1]*decision_variables[u2+i])-M_t[i][1]*decision_variables[a+i]-C_t[i][1]*(decision_variables[i]+decision_variables[i+1])/2
        return remainder
    return constraint1

def create_constraint2(n_discretization):
    a = n_discretization
    def constraint2(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = decision_variables[i+1]-decision_variables[i]-2*decision_variables[a+i]*1/(n_discretization-1)
        return remainder
    return constraint2

def create_b_bounds(n_discretization):
    lb=[]
    ub=[]
    lb.extend([1E-6]*n_discretization)
    lb.extend([-np.inf]*3*(n_discretization-1))
    ub.extend([np.inf]*(4*n_discretization-3))
    bounds = scp.optimize.Bounds(lb,ub)
    return bounds

def create_constraint3(mu,mass,n_discretization):
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    def constraint3(decision_variables):
        remainder = np.zeros(n_discretization-1)
        for i in range(n_discretization-1):
            remainder[i] = -(decision_variables[u1+i]**2+decision_variables[u2+i]**2)**0.5+mu*mass*9.81
        return remainder
    return constraint3

def build_x0(R_t,C_t,n_discretization):
    x0 = np.ones(n_discretization)*0.0001
    x0 = np.append(x0,np.zeros(3*(n_discretization-1)))
    u1=n_discretization+n_discretization-1
    u2=n_discretization+n_discretization-1+n_discretization-1
    for i in range(n_discretization-1):
        u = 0.0001*np.dot(np.linalg.inv(R_t[i]),C_t[i])
        x0[u1+i]=u[0]
        x0[u2+i]=u[1]
    return x0

def optimization(R_t,M_t,C_t,A_t,n_discretization,xsi):
    x0 =  build_x0(R_t,C_t,n_discretization)#{"u": np.ones((n_discretization,2)), "a": np.ones(n_discretization-1), "b": np.ones(n_discretization-1) }
    objective_function = create_objective(xsi, A_t,n_discretization)
    constraint1 = create_constraint1(R_t,M_t,C_t,n_discretization)
    constraint2=create_constraint2(n_discretization)
    mu=1
    mass=30
    constraint3 =create_constraint3(mu,mass,n_discretization)
    bounds = create_b_bounds(n_discretization)
    cons = [
    {'type': 'eq', 'fun': constraint1},  # Equality constraint 1
    {'type': 'eq', 'fun': constraint2},  # Equality constraint 2
    {'type': 'ineq', 'fun': constraint3},  # Inequality friction circle
        ]
    options = {
    'disp': True,      # Display iteration info
    'maxiter': 1000,   # Increase the maximum number of iterations
    'ftol': 1e-6,      # Tolerance on function value changes
        }       
    result = scp.optimize.minimize(objective_function, x0, method='SLSQP', constraints=cons,bounds=bounds,options=options)
    #print("constraint1",constraint1(result.x))
    #print("constraint2",constraint2(result.x))
    return  result, x0