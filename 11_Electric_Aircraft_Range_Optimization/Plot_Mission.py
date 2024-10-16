# Plot_Mission.py
# 
# Created: Mar, 2023. M. Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    
from RCAIDE.Library.Plots  import *    
# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(nexus): 
    results   = nexus.results.base
    
    # Plot Flight Conditions 
    plot_flight_conditions(results) 
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
    
    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results)
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results) 
    plot_battery_cell_conditions(results)
    plot_battery_degradation(results)
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_motor_and_rotor_efficiencies(results)
     
    return