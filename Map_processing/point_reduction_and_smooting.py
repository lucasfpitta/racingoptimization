import numpy as np
import scipy as sci






#project the point p1 from the smaller curve (less points) in to the
#longer curve array
def projection(p1,arr2):
    
    
    #find the closest point on the other array (only takes into 
    # account xy [0:2] coordinate)
    distance1,index1 = sci.spatial.KDTree(arr2[:,0:2]).query(p1[0:2])
    p_line2_1=arr2[index1]
    
    
    #finds the second closest point (only takes into 
    # account xy [0:2] coordinate)
    arr2 = np.delete(arr2,index1,0)
    distance2,index2 = sci.spatial.KDTree(arr2[:,0:2]).query(p1[0:2])
    p_line2_2=arr2[index2]
    
    
    #we project onto the line between the two closest points 
    if np.all(p_line2_1-p_line2_2) == 0:
        return p_line2_1
    
    #weight of the closest point
    t = np.dot(p1-p_line2_1,p_line2_2-p_line2_1)/np.dot(p_line2_2-p_line2_1
                                                    ,p_line2_2-p_line2_1)
    
    if t<0:
        new_point=p_line2_1
    elif t>1: 
        new_point=p_line2_2
    else:
        new_point = p_line2_1+t*(p_line2_2-p_line2_1)
    return(new_point)














#smoothing of the height error with mooving weighted average
def smooth_moving_average(array1,array2,n,alfa,max_slope):
    
    
    #creates the smoothed arrays
    smoothed_1, smoothed_2 = np.zeros(len(array1[0])), np.zeros(len(array2[0]))     
      
    
    for i in range(len(array1[0])):
        
         #if the points are the same, they're height is the same for continuity
        if array1[0,i] == array1[0,i-1] and array1[1,i] == array1[1,i-1]:
            smoothed_1[i] = smoothed_1[i-1]
        
        else:
            
            #index of neighborhood
            indices = range(i-int(n/2),i-int(n/2)+n) 
            
            #neighborhood points in both curves
            neighborhood_1 = array1[2].take(indices, mode = 'wrap')
            neighborhood_2 = array2[2].take(indices, mode = 'wrap')
            weights1, weights2 = np.zeros(n), np.zeros(n)
            
            #iterate over neighbors
            for j in range(n): 
                
                #point index (wrapping)
                k = i-int(n/2)+j-(i-int(n/2)+j)//len(array1[0])*len(array1[0]) 
                
                #defines the weight for each point
                #the weights are an exponencial decay over the xy cartesian distance
                # (alfa is the exp parameter)
                weights1[j] = np.exp(-alfa*(np.sqrt((array1[0,i]-array1[0,k]
                            )**2+(array1[1,i]-array1[1,k])**2+(array1[2,i]-
                                                        array1[2,k])**2))) 
                
                weights2[j] = np.exp(-alfa*(np.sqrt((array1[0,i]-array2[0,k]
                            )**2+(array1[1,i]-array2[1,k])**2+(array1[2,i]-
                                                        array2[2,k])**2)))
                
                
            #creates the smoothed array 1 with weighted average
            smoothed_1[i] = (np.dot(neighborhood_1, weights1)+np.dot(
                neighborhood_2, weights2))/np.sum(np.array([weights1,weights2])) 
            
            
            #height difference cap by the maximun slope of the track
            if np.abs(smoothed_1[i]-smoothed_1[i-1])/(np.sqrt((array1[0,i]-
                array1[0,i-1])**2+(array1[1,i]-array1[1,i-1])**2)
                                        ) > max_slope and i != 1: 
                smoothed_1[i] = smoothed_1[i-1]+np.sign((smoothed_1[i]-
                                smoothed_1[i-1]))*max_slope*np.sqrt((array1[0,i]-
                                array1[0,i-1])**2+(array1[1,i]-array1[1,i-1])**2)      
                                
                                
    smoothed_1[0]= smoothed_1[-1] #for continuity
    
    
    #repeat for the other side of the track
    for i in range(len(array2[0])):
        if array2[0,i] == array2[0,i-1] and array2[1,i] == array2[1,i-1]:
            smoothed_2[i] = smoothed_2[i-1]
            
            
        else:
            indices = range(i-int(n/2),i-int(n/2)+n)
            
            neighborhood_1 = array1[2].take(indices, mode = 'wrap')
            neighborhood_2 = array2[2].take(indices, mode = 'wrap')
            
            weights1, weights2 = np.zeros(n), np.zeros(n)
            
            for j in range(n):
                k = i-int(n/2)+j-(i-int(n/2)+j)//len(array2[0])*len(array2[0])
                weights1[j] = np.exp(-alfa*(np.sqrt((array2[0,i]-array1[0,k]
                            )**2+(array2[1,i]-array1[1,k])**2+(array2[2,i]-
                                                            array1[2,k])**2)))
                weights2[j] = np.exp(-alfa*(np.sqrt((array2[0,i]-array2[0,k]
                            )**2+(array2[1,i]-array2[2,k])**2+(array2[2,i]-
                                                            array2[2,k])**2)))
                
                
            smoothed_2[i] = (np.dot(neighborhood_1, weights1)+np.dot(
                neighborhood_2, weights2))/np.sum(np.array([weights1,weights2]))
            
            
            if np.abs(smoothed_2[i]-smoothed_2[i-1])/(np.sqrt((
                array2[0,i]-array2[0,i-1])**2+(array2[1,i]-array2[1,i-1]
                                        )**2)) > max_slope and i != 1:
                
                
                smoothed_2[i] = smoothed_2[i-1]+np.sign((smoothed_2[i]-
                                smoothed_2[i-1]))*max_slope*np.sqrt((
                                array2[0,i]-array2[0,i-1])**2+(array2[1,i]-
                                                        array2[1,i-1])**2)
                                
                                
    smoothed_2[0]= smoothed_2[-1]
    
    print("Track total elevation change", np.max(smoothed_1)-np.min(smoothed_1))
    
    return smoothed_1, smoothed_2
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#points reduction algorithm such that both arrays have the same length
def reduce_points(arr1,arr2):
    
    #arrays final length
    n = min(len(arr1[0]),len(arr2[0])) 
    
    #store all coordinates for both sides
    coords = np.zeros((6,n)) 
    
    
    #if arr1 is the smallest
    if len(arr1[0])<=len(arr2[0]):
        coords[0:3]=arr1
        
        #project every point
        for i in range(n): 
            [coords[3][i], coords[4][i], coords[5][i]]= projection(
                arr1[:,i],np.transpose(arr2))
            
    #if arr2 is the smallest        
    else:
        coords[3:6]=arr2
        
        #project every point
        for i in range(n): 
            [coords[0][i], coords[1][i], coords[2][i]] = projection(
                arr2[:,i],np.transpose(arr1))
            
            
    return(coords)
