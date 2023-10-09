'''

The script below documents how to set up and plot the results of a polar analysis of an airfoil 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
from RCAIDE.Core import Units
from RCAIDE.Methods.Aerodynamics.Airfoil_Panel_Method   import airfoil_analysis   
from RCAIDE.Methods.Geometry.Two_Dimensional.Airfoil    import  compute_naca_4series
from RCAIDE.Visualization  import *    
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():      
    # -----------------------------------------------
    # Batch analysis of single airfoil - NACA 2410 
    # -----------------------------------------------
    # angle of attack sweep at Reynolds number = 10,000
    
    # define Reynolds numbers 
    Re_vals            = np.atleast_2d(np.ones(7)*1E5)  
    
    # define angle of attack swee
    AoA_vals           = np.atleast_2d(np.linspace(-4,8,7)*Units.degrees) 
    
    # define airfoil name 
    airfoil_file       = '2410'
    
    # we are going to use the built-in geometry generator 
    airfoil_geometry   = compute_naca_4series(airfoil_file,npoints = 50)
    
    # compute airfoil properties using panel code 
    airfoil_properties = airfoil_analysis(airfoil_geometry,AoA_vals,Re_vals)  
    
    # Plots !     
    plot_airfoil_surface_forces(airfoil_properties)   
    plot_airfoil_polars(airfoil_properties) 
    plot_airfoil_boundary_layer_properties(airfoil_properties,show_legend = True)   
 
    return    

if __name__ == '__main__': 
    main()  
    plt.show()