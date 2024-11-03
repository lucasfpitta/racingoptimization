import numpy as np
import geopandas as gpd
import fiona
fiona.drvsupport.supported_drivers['kml'] = 'rw'
fiona.drvsupport.supported_drivers['KML'] = 'rw'



#function to read google_earth outlines and output outline in meters
#Inputs .kml documents
#Outputs 2d vectors of the outlines coordid_matrix 1 and coord_matrix2
def get_outline(doc1, doc2):

    #document reader
    #gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
    my_map1 = gpd.read_file(doc1, driver='KML')
    my_map2 = gpd.read_file(doc2, driver='KML')


    #extracting points
    my_map1['points'] = my_map1.apply(lambda x: [y for y in 
                            x['geometry'].coords], axis=1)
    my_map2['points'] = my_map2.apply(lambda x: [y for y in 
                            x['geometry'].coords], axis=1)


    #creates the coordinates matrix for each line
    coord_matrix1 = np.zeros((3,int(len(my_map1.points[0])+1)))
    coord_matrix2 = np.zeros((3,int(len(my_map2.points[0])+1)))
    
    
    #earth shape constants
    a=6378137
    b=6356752.3142
    e2=(a**2-b**2)/a**2
    
    
    
    #translate latitute and longitude to meters, with coordinated in 
    #respect to the first point. 

    for i in range(len(coord_matrix1[0])-1):
        coord_matrix1[0,i]=(my_map1.points[0][i][0]-my_map1.points[0][0][0]
            )*np.pi*a*np.cos(my_map1.points[0][0][0]*np.pi/180)/(180*np.sqrt(
            1-e2*np.sin(my_map1.points[0][0][0]*np.pi/180)))
        coord_matrix1[1,i]=(my_map1.points[0][i][1]-my_map1.points[0][0][1]
            )*np.pi*a*(1-e2)/(180*pow(1-e2*pow(np.sin(
            my_map1.points[0][0][1]*np.pi/180),2),3/2))
        coord_matrix1[2,i]=my_map1.points[0][i][2]-my_map1.points[0][0][2]
        
    #adds the last point which is the same as the first 
    coord_matrix1[0,-1]=coord_matrix1[0][0]
    coord_matrix1[1,-1]=coord_matrix1[1][0]
    coord_matrix1[2,-1]=coord_matrix1[2][0]
    
    
    #translate latitute and longitude to meters, with coordinated in 
    #respect to the first point. 
    
    for i in range(len(coord_matrix2[0])-1):
        coord_matrix2[0,i]=(my_map2.points[0][i][0]-my_map1.points[0][0][0]
            )*np.pi*a*np.cos(my_map1.points[0][0][0]*np.pi/180)/(180*np.sqrt(
            1-e2*np.sin(my_map1.points[0][0][0]*np.pi/180)))
        coord_matrix2[1,i]=(my_map2.points[0][i][1]-my_map1.points[0][0][1]
            )*np.pi*a*(1-e2)/(180*pow(1-e2*pow(np.sin(
            my_map1.points[0][0][1]*np.pi/180),2),3/2))
        coord_matrix2[2,i]=my_map2.points[0][i][2]-my_map1.points[0][0][2]
        
    #adds the last point which is the same as the first     
    coord_matrix2[0,-1]=coord_matrix2[0][0]
    coord_matrix2[1,-1]=coord_matrix2[1][0]
    coord_matrix2[2,-1]=coord_matrix2[2][0]
    
    
    
    
    return coord_matrix1, coord_matrix2





if __name__ == "__main__":
    right, left,  = get_outline('Map_processing/Maps_kml/extSEM.kml',
                                'Map_processing/Maps_kml/intSEM.kml')
    print(right,left)