# QS_tutorial.py
#  
#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units 
from RCAIDE.Library.Methods.Energy.Sources.Battery.Common             import initialize_from_mass 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion   import nasa_motor 
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor               import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor            import design_motor 
from RCAIDE.Library.Plots  import *    

import matplotlib.pyplot as plt
from copy import deepcopy
import os  

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    # vehicle data
    vehicle  = vehicle_setup() 

    # Set up vehicle configs
    configs  = configs_setup(vehicle)

    # create analyses
    analyses = analyses_setup(configs)

    # mission analyses 
    mission = mission_setup(analyses)

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
    results = missions.base_mission.evaluate()  

    # plot the results
    plot_mission(results)      
    
    return  

# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():
     #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------     

    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'tail_sitter'
 
    vehicle.mass_properties.takeoff         = 0.82 * Units.kg
    vehicle.mass_properties.operating_empty = 0.82 * Units.kg
    vehicle.mass_properties.max_takeoff     = 0.82 * Units.kg
    
    # basic parameters
    vehicle.reference_area                  = 0.1668 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------   


    wing = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.areas.reference         = vehicle.reference_area
    wing.spans.projected         = 1.03 * Units.m
    wing.aspect_ratio            = (wing.spans.projected**2)/wing.areas.reference 
    wing.sweeps.quarter_chord    = 5.0 * Units.deg
    wing.thickness_to_chord      = 0.12
    wing.taper                   = 1.0
    wing.dynamic_pressure_ratio  = 1.0
    wing.chords.mean_aerodynamic = 0.162 * Units.m
    wing.chords.root             = 0.162 * Units.m
    wing.chords.tip              = 0.162 * Units.m
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees
    wing.high_lift               = False
    wing.vertical                = False
    wing.symmetric               = True

    # add to vehicle
    vehicle.append_component(wing)  
    
  #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################## Energy Network ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    net                              = RCAIDE.Framework.Networks.Electric()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                              = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()
    bus.fixed_voltage                = False 
  
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Battery
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat                          = RCAIDE.Library.Components.Energy.Sources.Batteries.Lithium_Ion_Generic()
    bat.mass_properties.mass     = 0.17 * Units.kg
    bat.specific_energy          = 175.*Units.Wh/Units.kg
    bat.resistance               = 0.003
    bat.pack.nominal_voltage     = 11.1
    initialize_from_mass(bat)
    bus.batteries.append(bat)     
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Electronic Speed Controller    
    #------------------------------------------------------------------------------------------------------------------------------------  
    esc_1            = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc_1.tag        = 'esc_1'
    esc_1.efficiency = 0.95 
    bus.electronic_speed_controllers.append(esc_1)  
 
    esc_2            = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc_2.tag        = 'esc_2'
    esc_2.efficiency = 0.95 
    bus.electronic_speed_controllers.append(esc_2)    

    esc_3            = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc_3.tag        = 'esc_3'
    esc_3.efficiency = 0.95 
    bus.electronic_speed_controllers.append(esc_3)   

    esc_4            = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc_4.tag        = 'esc_4'
    esc_4.efficiency = 0.95 
    bus.electronic_speed_controllers.append(esc_4)   
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propeller    
    #------------------------------------------------------------------------------------------------------------------------------------           
    propeller                                   = RCAIDE.Library.Components.Propulsors.Converters.Propeller()  
    propeller.number_of_blades                  = 2.0
    propeller.tip_radius                        = 4.    * Units.inch
    propeller.hub_radius                        = 0.125 * Units.inch
    propeller.cruise.design_freestream_velocity = 15.0 # freestream m/s 
    propeller.cruise.design_angular_velocity    = 7500. * Units['rpm']
    propeller.cruise.design_Cl                  = 0.7
    propeller.cruise.design_altitude            = 0.1 * Units.km
    propeller.cruise.design_power               = 200. * Units.watts
    propeller                                   = design_propeller(propeller)   
    
    origins = [[0., 0.15, -0.05], [0., -0.15, -0.05], [0., .35, 0.05], [0., 0.35, 0.05]] 
    for ii in range(4):
        prop          = deepcopy(propeller)
        prop.tag      = 'propeller_'+ str(ii+1)
        prop.origin   = [origins[ii]] 
        bus.rotors.append(prop)  
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Motor
    #------------------------------------------------------------------------------------------------------------------------------------           
    motor                         = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor() 
    motor.efficiency              = 0.9
    motor.nominal_voltage         = bat.pack.nominal_voltage 
    motor.origin                  = propeller.origin 
    motor.propeller_radius        = propeller.tip_radius 
    motor.no_load_current         = 0.01  
    motor.wing_mounted            = True 
    motor.wing_tag                = 'main_wing'
    motor.rotor_radius            = propeller.tip_radius
    motor.design_torque           = propeller.cruise.design_torque
    motor.angular_velocity        = propeller.cruise.design_angular_velocity/motor.gear_ratio  
    motor                         = design_motor(motor)
    motor.mass_properties.mass    = nasa_motor(motor.design_torque)    
    

    for ii in range(4):
        rotor_motor = deepcopy(motor)
        rotor_motor.tag    = 'motor_' + str(ii+1)
        rotor_motor.origin = [origins[ii]]
        bus.motors.append(rotor_motor)  
    
    # append bus   
    net.busses.append(bus)        
    
    # Component 4 the Payload
    payload = RCAIDE.Library.Components.Payloads.Payload()
    payload.power_draw           = 0. #Watts 
    payload.mass_properties.mass = 0.0 * Units.kg
    net.payload                  = payload
    
    # Component 5 the Avionics
    avionics                     = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw          = 2. #Watts  
    net.avionics                 = avionics      
 
    vehicle.append_energy_network(net)  

    return vehicle


# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)
    
    return configs

 

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
    
    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Framework.Analyses.Weights.Weights()
    weights.settings.empty_weight_method =  RCAIDE.Library.Methods.Weights.Correlation_Buildups.UAV.empty
    weights.vehicle = vehicle
    analyses.append(weights)
    
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    aerodynamics.settings.maximum_lift_coefficient   = 1.5
    analyses.append(aerodynamics)    
    
    # ------------------------------------------------------------------
    #  Energy
    energy = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle  = vehicle 
    analyses.append(energy)
    
    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)
    
    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   
    
    # done!
    return analyses    

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------
def mission_setup(analyses):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Framework.Mission.Sequential_Segments()  
    mission.atmosphere  = RCAIDE.Library.Attributes.Atmospheres.Earth.US_Standard_1976()
    mission.planet      = RCAIDE.Library.Attributes.Planets.Earth()
    
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments
    
    # base segment
    base_segment = Segments.Segment()    
    base_segment.state.numerics.number_control_points  = 3

    #------------------------------------------------------------------    
    #  Climb Hover
    #------------------------------------------------------------------    
    
    segment = RCAIDE.Framework.Mission.Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "Climb" 
    segment.analyses.extend(analyses.base) 
    ones_row                                   = segment.state.ones_row 
    segment.initial_battery_state_of_charge    = 0.89 
    segment.altitude_start                     = 0.
    segment.altitude_end                       = 100. * Units.m
    segment.climb_rate                         = 3.  * Units.m / Units.s 
    segment.air_speed                          = 3.  * Units.m / Units.s
    segment.state.unknowns.body_angle          = ones_row(1) * 90. *Units.deg  
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------    
    #   Hover
    # ------------------------------------------------------------------     
    segment                                                      = RCAIDE.Framework.Mission.Segments.Vertical_Flight.Hover(base_segment)
    segment.tag                                                  = "Hover_1" 
    segment.analyses.extend(analyses.base) 
    segment.time                                                 = 60* Units.seconds
    segment.state.conditions.frames.body.inertial_rotations[:,1] = ones_row(1) * 90.*Units.deg    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------    
    #   Hover Transition
    # ------------------------------------------------------------------   
    segment                   = RCAIDE.Framework.Mission.Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag               = "Transition_to_Cruise" 
    segment.analyses.extend(analyses.base) 
    segment.initial_battery_state_of_charge    = 0.89 
    segment.acceleration      = 1.5 * Units['m/s/s']
    segment.air_speed_start   = 0.0
    segment.air_speed_end     = 15.0 
    segment.altitude          = 100. * Units.m  
    mission.append_segment(segment)   

    # ------------------------------------------------------------------    
    #   Cruise
    # ------------------------------------------------------------------      
    segment                                 = RCAIDE.Framework.Mission.Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                             = "Cruise" 
    segment.analyses.extend(analyses.base) 
    segment.distance                        = 3.  * Units.km
    segment.air_speed                       = 15. * Units.m/Units.s
    segment.altitude                        = 100. * Units.m   
    mission.append_segment(segment)            
    
    # ------------------------------------------------------------------    
    #   Hover Transition
    # ------------------------------------------------------------------     
    
    segment                                 = RCAIDE.Framework.Mission.Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                             = "Transition_to_hover" 
    segment.analyses.extend(analyses.base) 
    segment.acceleration                    = -0.5 * Units['m/s/s'] 
    segment.air_speed_end                   = 0.0 
    segment.altitude                        = 100. * Units.m  
    
    mission.append_segment(segment)  
    
    # ------------------------------------------------------------------    
    #   Hover
    # ------------------------------------------------------------------     
    segment                                                      = RCAIDE.Framework.Mission.Segments.Vertical_Flight.Hover(base_segment)
    segment.tag                                                  = "Hover_2" 
    segment.analyses.extend(analyses.base)                      
    segment.time                                                 = 60* Units.seconds
    segment.state.conditions.frames.body.inertial_rotations[:,1] = ones_row(1)* 90.*Units.deg  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------    
    #   Descent Hover
    # ------------------------------------------------------------------      
    segment              = RCAIDE.Framework.Mission.Segments.Vertical_Flight.Descent(base_segment)
    segment.tag          = "Descent" 
    segment.analyses.extend(analyses.base) 
    segment.altitude_end = 0. * Units.m
    segment.descent_rate = 3. * Units.m / Units.s   
    segment.state.conditions.frames.body.inertial_rotations[:,1] = ones_row(1)* 90.*Units.deg  
    mission.append_segment(segment)       
    
    #------------------------------------------------------------------    
    #   Mission definition complete    
    #-------------------------------------------------------------------
    
    return mission


def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def missions_setup(mission): 

    missions         = RCAIDE.Framework.Mission.Missions()

    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions  

# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------
def plot_mission(results):


    # Plot Flight Conditions 
    plot_flight_conditions(results) 
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
    
    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results)
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results)

    # Plot propeller Disc and Power Loading
    plot_rotor_conditions(results)  
    
    return     

if __name__ == '__main__':
    main()
    plt.show()