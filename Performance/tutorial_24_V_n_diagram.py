# Regression/scripts/Tests/performance_payload_range.py
# 
# 
# Created:  Jul 2023, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports  
import RCAIDE
from RCAIDE.Framework.Core import Units , Container
from RCAIDE.Library.Methods.Performance.compute_payload_range_diagram        import compute_payload_range_diagram

# python imports     
import numpy as np  
import sys
import matplotlib.pyplot as plt  
import os
# local imports 
sys.path.append(os.path.join( os.path.split(os.path.split(sys.path[0])[0])[0], 'Vehicles')) 
# ----------------------------------------------------------------------------------------------------------------------
#  REGRESSION
# ----------------------------------------------------------------------------------------------------------------------  
def main():
    
    # vehicle data
    vehicle             = vehicle_setup()

    # take out control surfaces to make regression run faster
    for wing in vehicle.wings:
        wing.control_surfaces  = Container()
        
    assigned_propulsors = [['starboard_propulsor','port_propulsor']]   
    altitude            = 15000   * Units.feet 
    airspeed            = 130 * Units.kts
    max_range_guess     =  20.   * Units.nautical_mile
    

    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle 

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics          = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle  = vehicle
    aerodynamics.settings.number_of_spanwise_vortices   = 5
    aerodynamics.settings.number_of_chordwise_vortices  = 2  
    
    # run payload range analysis 
    payload_range_results =  compute_payload_range_diagram(vehicle,assigned_propulsors,airspeed, altitude,max_range_guess)
                                   
    return payload_range_results  


if __name__ == '__main__': 
    main()    
    plt.show() 