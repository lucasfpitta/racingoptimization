import numpy as np


#Calculates the time to traverse each section using the generalized velocity squared b
#Input 1d b vector
#Output 1d time vector
def reconstruct(b):
    dimension = len(b)
    t=np.zeros(dimension)
    for i in range(dimension-1):
        t[i+1] = 2*1/(dimension-1)/(b[i]**0.5+b[i+1]**0.5)+t[i]
    return t
