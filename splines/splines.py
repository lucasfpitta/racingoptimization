import numpy as np
import scipy as scp

def knot_points(left, right, alfas):
    splpoints = np.zeros((2,len(right[0])))
    for i in range(len(right[0])):
        splpoints[0][i] = right[0][i]+alfas[i]*(left[0][i]-right[0][i])
        splpoints[1][i] = right[1][i]+alfas[i]*(left[1][i]-right[1][i])
    return splpoints

def splines_and_derivatives(splpoints,N):
    t = np.linspace(0,1,num = len(splpoints[0]))
    interpol = scp.interpolate.CubicSpline(t, (splpoints[0],splpoints[1]),axis=1, bc_type='natural')
    #spline = interpol(np.linspace(0,1,num = N))
    diff = interpol.derivative()
    #derivative = diff(np.linspace(0,1,num = N))
    return interpol, diff#spline, derivative

def find_angle(derivative):
    angle = np.zeros(len(derivative[0]))
    for i in range(len(derivative[0])):
        if derivative[0][i] >= 0:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))
        else:
            angle[i] = np.arctan((derivative[1][i])/(derivative[0][i]))+np.pi
    angle[0]=angle[-1]
    return angle

def path_info(left, right, alfas,N):
    splpoints = knot_points(left, right, alfas)
    spline, derivative = splines_and_derivatives(splpoints,N)
    angle = find_angle(derivative(np.linspace(0,1,num = N)))
    return spline, derivative, angle