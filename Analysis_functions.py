##################################
####          PACKAGES        ####
##################################
import numpy as np
from scipy.interpolate import griddata

#################################
####        FUNCTIONS        ####
#################################

def get_volume(x, y, LB_idx, UB_idx):
    volume_trapz = np.abs(np.trapz(y[UB_idx:LB_idx]-y[LB_idx], x[UB_idx:LB_idx]))

    return volume_trapz

def get_gradient(x, y, years, seaward_bound, landward_bound):
    gradient = []
    for i, yr in enumerate(years):
        # Extract elevation profile with seaward and landward boundaries
        elevation_y = []
        cross_shore_x = []
        
        for xc in range(len(x)): 
            if x[xc] < seaward_bound[i] and x[xc] > landward_bound[i] and np.isnan(y[xc,i]) == False:
                elevation_y.append(y[xc,i])
                cross_shore_x.append(x[xc])
        # Calculate gradient for domain
        if cross_shore_x == []:
            gradient.append(np.nan)
        else:
            gradient.append(np.polyfit(cross_shore_x, elevation_y, 1)[0])
        
    return gradient

def find_intersections(x, y, y_value):
    value_vec = []
    for x_val in range(len(x)):
        value_vec.append(y_value)
    
    diff = np.nan_to_num(np.diff(np.sign(y - value_vec)))
    intersections = np.nonzero(diff)
    
    return intersections

