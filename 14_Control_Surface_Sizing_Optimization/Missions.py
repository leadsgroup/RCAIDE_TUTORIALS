# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import Units 
import numpy as np

#   Define the Mission
def stick_fixed_stability_setup(analyses,vehicle): 
    missions                     = RCAIDE.Framework.Mission.Mission.Container()  
    max_speed_multiplier         = 1.0 # this multiplier is used to compute V_max from V_nominal
    missions.stick_fixed_cruise  = base_mission_setup(vehicle,max_speed_multiplier) 
 
    return missions   

def elevator_sizing_setup(analyses,vehicle): 
    missions = RCAIDE.Framework.Mission.Mission.Container() 
    max_speed_multiplier      = 1.4 # this multiplier is used to compute V_max from V_nominal
    missions.elevator_sizing  = base_mission_setup(vehicle,max_speed_multiplier)   
 
    return missions   

def aileron_rudder_sizing_setup(analyses,vehicle): 
    missions = RCAIDE.Framework.Mission.Mission.Container() 
    max_speed_multiplier      = 1.0     
    missions.aileron_sizing   = base_mission_setup(vehicle,max_speed_multiplier)  
    max_speed_multiplier      = 1.4   # this multiplier is used to compute V_max from V_nominal   
    missions.turn_criteria    = base_mission_setup(vehicle,max_speed_multiplier) 
 
    return missions   
    
def flap_sizing_setup(analyses,vehicle): 
    missions = RCAIDE.Framework.Mission.Mission.Container() 
    max_speed_multiplier     = 1.0      
    missions.flap_sizing     = base_mission_setup(vehicle,max_speed_multiplier)   
    return missions        
    

# ------------------------------------------------------------------
#   Initialize the Mission
# ------------------------------------------------------------------    
    
def base_mission_setup(vehicle,max_speed_multiplier):   
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # airport
    airport = RCAIDE.Library.Attributes.Airports.Airport()
    airport.altitude   =  0. * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = RCAIDE.Library.Attributes.Atmospheres.Earth.US_Standard_1976()

    mission.airport = airport    
 
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    # base segment
    base_segment = Segments.Segment()
    ones_row     = base_segment.state.ones_row
    base_segment.state.numerics.number_control_points    = 3
 
    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                       = "cruise"   
    segment.altitude                  = 8012   * Units.feet
    segment.air_speed                 = 120.91 * Units['mph'] * max_speed_multiplier
    segment.distance                  =  20.   * Units.nautical_mile   
    segment.state.unknowns.throttle   = 0.5 *  ones_row(1)   
    mission.append_segment(segment)     
    
    return mission
