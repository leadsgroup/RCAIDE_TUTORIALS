# Missions.py
# 
# Created:  Mar 2023, M. Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------  
import RCAIDE
from RCAIDE.Core import Units     
from RCAIDE.Methods.Performance.estimate_stall_speed   import estimate_stall_speed 
 

# ----------------------------------------------------------------------
#   Missions Setup
# ---------------------------------------------------------------------- 
    
def setup(analyses):
    
    # the mission container
    missions = RCAIDE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------
    base_mission = mission_setup(analyses)
    missions.base = base_mission 
 
    return missions  
 

#------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def mission_setup(analyses): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission'

    # airport  
    airport            = RCAIDE.Attributes.Airports.Airport()   
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport      

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = airport.altitude,temperature_deviation= 1.)     
    mission.airport    = airport       

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                             = Segments.Segment() 
    ones_row                                                 = base_segment.state.ones_row 
    base_segment                                             = Segments.Segment()   
    base_segment.state.numerics.number_control_points        = 8
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)         

    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.1  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0   
    segment.battery_energy                                   =  analyses.base.energy.networks.all_electric.battery.pack.max_energy    
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min']  
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    mission.append_segment(segment)  
                
    #------------------------------------------------------------------
    #  Initial Climb Area Segment Flight 1  
    #------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']           
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])  
    mission.append_segment(segment)  
     
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'       
    segment.analyses.extend( analyses.base )         
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 2500 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = 175.* Units['mph']    
    segment.climb_rate                                       = 500* Units['ft/min']      
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)  
    
    vehicle        = analyses.base.aerodynamics.geometry 
    # ------------------------------------------------------------------
    #   Cruise Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Cruise'  
    segment.tag = 'cruise'  
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 175.*Units['mph']   
    segment.altitude                                         = 2500 * Units.feet 
    segment.air_speed                                        = 175.*Units['mph']   
    segment.altitude                                         = 2500 * Units.feet 
    segment.distance                                         = 50*Units.nmi           
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)   
 
    
    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Descent'   
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 2500 * Units.feet  
    segment.altitude_end                                     = 1000.0 * Units.feet
    segment.air_speed_start                                  = 175.* Units['mph']    
    segment.air_speed_end                                    = Vstall*1.3     
    segment.climb_rate                                       = -300 * Units['ft/min']    
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)  
     

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg' 
    segment.analyses.extend(analyses.base) 
    segment.air_speed_start                                  = Vstall*1.3 
    segment.air_speed_end                                    = Vstall*1.2       
    segment.distance                                         =  6000 * Units.feet
    segment.acceleration                                     = -0.05 * Units['m/s/s']     
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.7])  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg' 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.2
    segment.air_speed_end                                    = Vstall*1.1  
    segment.climb_rate                                       = -300 * Units['ft/min']  
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment ,  initial_throttles = [0.7])
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)    
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)    
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    segment.battery_energy                                   =  analyses.base.energy.networks.all_electric.battery.pack.max_energy    
    segment =  analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)   
    
    
    return mission 
