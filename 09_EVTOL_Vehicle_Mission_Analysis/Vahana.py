''' 
# Tiltwing_EVTOL.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
import RCAIDE
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common                    import initialize_from_circuit_configuration 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion            import compute_motor_weight
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor                     import design_motor
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor                        import design_prop_rotor ,design_prop_rotor 
from RCAIDE.Library.Methods.Weights.Physics_Based_Buildups.Electric            import converge_physics_based_weight_buildup 
from RCAIDE.Library.Plots                                                      import *      
# python imports 
import os
import numpy as np 
from copy import deepcopy
import pickle
import  pandas as pd
import matplotlib.pyplot as plt  

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main(): 
    # vehicle data 
    vehicle  = vehicle_setup() 
        
    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate() 
     
    # plot the results 
    plot_results(results)
      
    return 

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------
def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses


def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
    
    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics          = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)   

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
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
def vehicle_setup(new_regression=True): 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    vehicle                                     = RCAIDE.Vehicle()
    vehicle.tag                                 = 'Vahana'
    vehicle.configuration                       = 'eVTOL'
         
    # mass properties
    vehicle.mass_properties.takeoff             = 735. 
    vehicle.mass_properties.operating_empty     = 735.
    vehicle.mass_properties.max_takeoff         = 735.
    vehicle.mass_properties.center_of_gravity   = [[ 2.0144,   0.  ,  0.]] 
    vehicle.passengers                          = 0
    vehicle.envelope.ultimate_load              = 5.7
    vehicle.envelope.limit_load                 = 3.     

    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------
    wing                                        = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                                    = 'canard_wing'  
    wing.aspect_ratio                           = 11.37706641  
    wing.sweeps.quarter_chord                   = 0.0
    wing.thickness_to_chord                     = 0.18  
    wing.taper                                  = 1.  
    wing.spans.projected                        = 6.65 
    wing.chords.root                            = 0.95 
    wing.total_length                           = 0.95   
    wing.chords.tip                             = 0.95 
    wing.chords.mean_aerodynamic                = 0.95   
    wing.dihedral                               = 0.0  
    wing.areas.reference                        = wing.chords.root*wing.spans.projected 
    wing.areas.wetted                           = 2*wing.chords.root*wing.spans.projected*0.95  
    wing.areas.exposed                          = 2*wing.chords.root*wing.spans.projected*0.95 
    wing.twists.root                            = 0.  
    wing.twists.tip                             = 0.  
    wing.origin                                 = [[0.1,  0.0 , 0.0]]  
    wing.aerodynamic_center                     = [0., 0., 0.]     
    wing.winglet_fraction                       = 0.0 
    wing.symmetric                              = True

    
    ospath                                      = os.path.abspath(__file__) 
    separator                                   = os.path.sep
    rel_path                                    = os.path.dirname(ospath) + separator + '..'    + separator 
    airfoil                                     = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.coordinate_file                     = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'
    
    wing.append_airfoil(airfoil)
                                                
    # add to vehicle                                          
    vehicle.append_component(wing)                            
                                                
    wing                                        = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                                    = 'main_wing'  
    wing.aspect_ratio                           = 11.37706641  
    wing.sweeps.quarter_chord                   = 0.0
    wing.thickness_to_chord                     = 0.18  
    wing.taper                                  = 1.  
    wing.spans.projected                        = 6.65 
    wing.chords.root                            = 0.95 
    wing.total_length                           = 0.95   
    wing.chords.tip                             = 0.95 
    wing.chords.mean_aerodynamic                = 0.95   
    wing.dihedral                               = 0.0  
    wing.areas.reference                        = wing.chords.root*wing.spans.projected 
    wing.areas.wetted                           = 2*wing.chords.root*wing.spans.projected*0.95  
    wing.areas.exposed                          = 2*wing.chords.root*wing.spans.projected*0.95 
    wing.twists.root                            = 0.  
    wing.twists.tip                             = 0.  
    wing.origin                                 = [[ 5.138, 0.0  ,  1.323 ]]  # for images 1.54
    wing.aerodynamic_center                     = [0., 0., 0.]     
    wing.winglet_fraction                       = 0.0  
    wing.symmetric                              = True  
    vehicle.reference_area                      = 2*wing.areas.reference 
    wing.append_airfoil(airfoil)

    # add to vehicle 
    vehicle.append_component(wing)   


    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Fuselage ############################################################## 
    #------------------------------------------------------------------------------------------------------------------------------------
    
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage' 
    fuselage.seats_abreast                      = 0.  
    fuselage.seat_pitch                         = 1.  
    fuselage.fineness.nose                      = 1.5 
    fuselage.fineness.tail                      = 4.0 
    fuselage.lengths.nose                       = 1.7   
    fuselage.lengths.tail                       = 2.7 
    fuselage.lengths.cabin                      = 1.7  
    fuselage.lengths.total                      = 6.1  
    fuselage.width                              = 1.15  
    fuselage.heights.maximum                    = 1.7 
    fuselage.heights.at_quarter_length          = 1.2  
    fuselage.heights.at_wing_root_quarter_chord = 1.7  
    fuselage.heights.at_three_quarters_length   = 0.75 
    fuselage.areas.wetted                       = 12.97989862  
    fuselage.areas.front_projected              = 1.365211404  
    fuselage.effective_diameter                 = 1.318423736  
    fuselage.differential_pressure              = 0.  

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_0'   
    segment.percent_x_location                  = 0.  
    segment.percent_z_location                  = 0.  
    segment.height                              = 0.09  
    segment.width                               = 0.23473  
    segment.length                              = 0.  
    segment.effective_diameter                  = 0. 
    fuselage.Segments.append(segment)             

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'   
    segment.percent_x_location                  = 0.97675/6.1 
    segment.percent_z_location                  = 0.21977/6.1
    segment.height                              = 0.9027  
    segment.width                               = 1.01709  
    fuselage.Segments.append(segment)             


    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'    
    segment.percent_x_location                  = 1.93556/6.1 
    segment.percent_z_location                  = 0.39371/6.1
    segment.height                              = 1.30558   
    segment.width                               = 1.38871  
    fuselage.Segments.append(segment)             


    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'    
    segment.percent_x_location                  = 3.44137/6.1 
    segment.percent_z_location                  = 0.57143/6.1
    segment.height                              = 1.52588 
    segment.width                               = 1.47074 
    fuselage.Segments.append(segment)             

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 4.61031/6.1
    segment.percent_z_location                  = 0.10893
    segment.height                              = 1.3906
    segment.width                               = 1.11463  
    fuselage.Segments.append(segment)              

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 0.9827
    segment.percent_z_location                  = 0.180
    segment.height                              = 0.6145
    segment.width                               = 0.3838
    fuselage.Segments.append(segment)            
    
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 1. 
    segment.percent_z_location                  = 0.2058
    segment.height                              = 0.4
    segment.width                               = 0.25
    fuselage.Segments.append(segment)        

    # add to vehicle
    vehicle.append_component(fuselage)    
   
    sys                            = RCAIDE.Library.Components.Systems.System()
    sys.mass_properties.mass       = 5 # kg   
    vehicle.append_component(sys)    

    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################  Energy Network  ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------
    # define network
    network                                                = RCAIDE.Framework.Networks.Electric() 
    network.charging_power                                 = 1000
    #==================================================================================================================================== 
    # Lift Bus 
    #====================================================================================================================================          
    bus                                                    = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()
    bus.tag                                                = 'bus' 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus Battery
    #------------------------------------------------------------------------------------------------------------------------------------ 
    bat                                                    = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_NMC()
    number_of_modules                                      = 10 
    bat.tag                                                = 'bus_battery'
    bat.electrical_configuration.series                     = 8 
    bat.electrical_configuration.parallel                   = 60
    bat.cell.maximum_voltage                                = 4.2                                                                          
    bat.cell.nominal_capacity                               = 3.0                                                                          
    bat.cell.nominal_voltage                                = 3.6                                                                          
    initialize_from_circuit_configuration(bat)  
   
    bat.geometrtic_configuration.total                      = bat.electrical_configuration.total
    bat.voltage                                             = bat.maximum_voltage 
    bat.geometrtic_configuration.normal_count               = 20
    bat.geometrtic_configuration.parallel_count             = 24  
    
    for _ in range(number_of_modules):
        bus.battery_modules.append(deepcopy(bat))  
    
    for battery_module in  bus.battery_modules:
        bus.voltage  +=   battery_module.voltage 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Lift Propulsors 
    #------------------------------------------------------------------------------------------------------------------------------------    
     
    # Define Lift Propulsor Container 
    lift_propulsor                                = RCAIDE.Library.Components.Propulsors.Electric_Rotor()
    lift_propulsor.tag                            = 'lift_propulsor'      
    lift_propulsor.wing_mounted                   = True 
              
    # Electronic Speed Controller           
    prop_rotor_esc                                = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    prop_rotor_esc.efficiency                     = 0.95    
    prop_rotor_esc.tag                            = 'prop_rotor_esc_1'  
    lift_propulsor.electronic_speed_controller    = prop_rotor_esc  
    
    # Lift Rotor Design
    g                                             = 9.81                                    # gravitational acceleration   
    Hover_Load                                    = vehicle.mass_properties.takeoff*g *1.1  # hover load   

    prop_rotor                                    = RCAIDE.Library.Components.Propulsors.Converters.Prop_Rotor()   
    prop_rotor.tag                                = 'prop_rotor'   
    prop_rotor.tip_radius                         = 0.8875
    prop_rotor.hub_radius                         = 0.10 * prop_rotor.tip_radius
    prop_rotor.number_of_blades                   = 3
    prop_rotor.hover.design_altitude              = 40 * Units.feet   
    prop_rotor.hover.design_thrust                = Hover_Load/8 
    prop_rotor.hover.design_freestream_velocity   = np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))
    prop_rotor.hover.design_angular_velocity      = 0.65 * 343 /prop_rotor.tip_radius
    prop_rotor.hover.design_power_coefficient     = 0.02 # Guess 
    prop_rotor.oei.design_altitude                = 40 * Units.feet  
    prop_rotor.oei.design_thrust                  = Hover_Load/7  
    prop_rotor.oei.design_freestream_velocity     = np.sqrt(prop_rotor.oei.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))   
    prop_rotor.cruise.design_altitude             = 1500 * Units.feet
    prop_rotor.cruise.design_thrust               = 200    
    prop_rotor.cruise.design_freestream_velocity  = 130.* Units['mph']  
    
    airfoil                                       = RCAIDE.Library.Components.Airfoils.Airfoil()   
    airfoil.coordinate_file                       =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                           = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                     rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    prop_rotor.append_airfoil(airfoil)                
    prop_rotor.airfoil_polar_stations             = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]    
    prop_rotor.fidelity                          = "Momentum_Theory"
    design_prop_rotor(prop_rotor)  
    lift_propulsor.rotor =  prop_rotor
    
    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Motor  
    #------------------------------------------------------------------------------------------------------------------------------------    
    prop_rotor_motor                         = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    prop_rotor_motor.efficiency              = 0.95
    prop_rotor_motor.nominal_voltage         = bus.voltage * 0.75
    prop_rotor_motor.prop_rotor_radius       = prop_rotor.tip_radius 
    prop_rotor_motor.no_load_current         = 0.1  
    prop_rotor_motor.rotor_radius            = prop_rotor.tip_radius
    prop_rotor_motor.design_torque           = (prop_rotor.hover.design_thrust * np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))) /prop_rotor.hover.design_angular_velocity 
    prop_rotor_motor.angular_velocity        = prop_rotor.hover.design_angular_velocity/prop_rotor_motor.gear_ratio  
    design_motor(prop_rotor_motor)
    prop_rotor_motor.speed_constant          = 1.119848656457751 
    prop_rotor_motor.resistance              = 0.12534232206970
    prop_rotor_motor.mass_properties.mass    = compute_motor_weight(prop_rotor_motor.design_torque)     
    lift_propulsor.motor                     = prop_rotor_motor
     

    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Nacelle
    #------------------------------------------------------------------------------------------------------------------------------------     
    nacelle                           = RCAIDE.Library.Components.Nacelles.Nacelle() 
    nacelle.length                    = 0.45
    nacelle.diameter                  = 0.3 
    nacelle.flow_through              = False    
    lift_propulsor.nacelle            =  nacelle  

    # Front Rotors Locations 
    origins = [[-0.2, 1.347, 0.0], [-0.2, 3.2969999999999997, 0.0], [-0.2, -1.347, 0.0], [-0.2, -3.2969999999999997, 0.0],\
               [4.938, 1.347, 1.54], [4.938, 3.2969999999999997, 1.54],[4.938, -1.347, 1.54], [4.938, -3.2969999999999997, 1.54]] 
    
    for i in range(8): 
        lift_propulsor_i                                       = deepcopy(lift_propulsor)
        lift_propulsor_i.tag                                   = 'lift_rotor_propulsor_' + str(i + 1)
        lift_propulsor_i.rotor.tag                             = 'prop_rotor_' + str(i + 1) 
        lift_propulsor_i.rotor.origin                          = [origins[i]]  
        lift_propulsor_i.motor.tag                             = 'prop_rotor_motor_' + str(i + 1)  
        if i < 4: 
            lift_propulsor_i.motor.wing_tag                = 'canard_wing'
        else:
            lift_propulsor_i.motor.wing_tag                = 'main_wing'
        lift_propulsor_i.motor.origin                          = [origins[i]]  
        lift_propulsor_i.electronic_speed_controller.tag       = 'prop_rotor_esc_' + str(i + 1)  
        lift_propulsor_i.electronic_speed_controller.origin    = [origins[i]]  
        lift_propulsor_i.nacelle.tag                           = 'prop_rotor_nacelle_' + str(i + 1)  
        lift_propulsor_i.nacelle.origin                        = [origins[i]]   
        bus.propulsors.append(lift_propulsor_i)  

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Additional Bus Loads
    #------------------------------------------------------------------------------------------------------------------------------------            
    # Payload   
    payload                         = RCAIDE.Library.Components.Systems.Avionics()
    payload.power_draw              = 10. # Watts 
    payload.mass_properties.mass    = 1.0 * Units.kg
    bus.payload                     = payload 
                             
    # Avionics                            
    avionics                        = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw             = 10. # Watts  
    avionics.mass_properties.mass   = 1.0 * Units.kg
    bus.avionics                    = avionics    
   
    network.busses.append(bus)
    
    # append motor origin spanwise locations onto wing data structure
    motor_origins_front                                   = np.array(origins[:4])
    motor_origins_rear                                    = np.array(origins[5:])
    vehicle.wings['canard_wing'].motor_spanwise_locations = motor_origins_front[:,1]/ vehicle.wings['canard_wing'].spans.projected
    vehicle.wings['canard_wing'].motor_spanwise_locations = motor_origins_front[:,1]/ vehicle.wings['canard_wing'].spans.projected
    vehicle.wings['main_wing'].motor_spanwise_locations   = motor_origins_rear[:,1]/ vehicle.wings['main_wing'].spans.projected
      
        
    # append energy network 
    vehicle.append_energy_network(network)  
 
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------   
    converged_vehicle, breakdown = converge_physics_based_weight_buildup(vehicle)  
    print(breakdown) 

    return converged_vehicle

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    '''
    The configration set up below the scheduling of the nacelle angle and vehicle speed.
    Since one prop_rotor operates at varying flight conditions, one must perscribe  the 
    pitch command of the prop_rotor which us used in the variable pitch model in the analyses
    Note: low pitch at take off & low speeds, high pitch at cruise
    '''
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config                                                       = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag                                                   = 'base'     
    configs.append(base_config) 
 
    # ------------------------------------------------------------------
    #   Hover Climb Configuration
    # ------------------------------------------------------------------
    config                                                 = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                             = 'vertical_climb'
    vector_angle                                           = 90.0 * Units.degrees
    config.wings.main_wing.twists.root                     = vector_angle
    config.wings.main_wing.twists.tip                      = vector_angle
    config.wings.canard_wing.twists.root                   = vector_angle
    config.wings.canard_wing.twists.tip                    = vector_angle    
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
    configs.append(config)

    # ------------------------------------------------------------------
    #    
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    vector_angle                                      = 30.0  * Units.degrees 
    config.tag                                        = 'vertical_transition'
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
                propulsor.rotor.pitch_command   = propulsor.rotor.hover.design_pitch_command * 0.5 
    configs.append(config) 

    # ------------------------------------------------------------------
    #   Hover-to-Cruise Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                        = 'climb_transition'
    vector_angle                                      = 5.0  * Units.degrees  
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle 
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
                propulsor.rotor.pitch_command     = propulsor.rotor.cruise.design_pitch_command  
    configs.append(config) 

    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                        = 'cruise'   
    vector_angle                                      = 0.0 * Units.degrees 
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle  
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
                propulsor.rotor.pitch_command   = propulsor.rotor.cruise.design_pitch_command  
    configs.append(config)     
    
    # ------------------------------------------------------------------
    #   
    # ------------------------------------------------------------------ 
    config                                                 = RCAIDE.Library.Components.Configs.Config(vehicle)
    vector_angle                                           = 75.0  * Units.degrees   
    config.tag                                             = 'descent_transition'   
    config.wings.main_wing.twists.root                     = vector_angle
    config.wings.main_wing.twists.tip                      = vector_angle
    config.wings.canard_wing.twists.root                   = vector_angle
    config.wings.canard_wing.twists.tip                    = vector_angle
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
                propulsor.rotor.pitch_command   = propulsor.rotor.cruise.design_pitch_command * 0.5
    configs.append(config)  


    # ------------------------------------------------------------------
    #   Approach Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                        = 'approach_transition'   
    vector_angle                                      = 0.0 * Units.degrees 
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle 
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0] 
                propulsor.rotor.pitch_command   = propulsor.rotor.cruise.design_pitch_command *0.5  
    configs.append(config)
    
    
    # ------------------------------------------------------------------
    #   Hover Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                        = 'vertical_descent'
    vector_angle                                      = 90.0  * Units.degrees   
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle     
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors:
                propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
    configs.append(config)

    return configs 
    

def mission_setup(analyses ):
    

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment() 
     
  
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                                          = Segments.Vertical_Flight.Hover(base_segment)
    segment.tag                                                      = "Hover"   
    segment.analyses.extend(analyses.vertical_climb)
    
    segment.altitude                                                 = 40.  * Units.ft  
    segment.initial_battery_state_of_charge                          = 1.0 
                        
    # define flight dynamics to model              
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                            'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
      
    mission.append_segment(segment)
 
 

    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                                          = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                                      = "Vertical_Climb_1"   
    segment.analyses.extend(analyses.vertical_climb)                
    segment.altitude_end                                             = 60.  * Units.ft  
    segment.initial_battery_state_of_charge                          = 1.0 
    segment.climb_rate                                               = 100. * Units['ft/min'] 
          
    # define flight dynamics to model            
    segment.flight_dynamics.force_z                                  = True 

    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                            'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    
    mission.append_segment(segment)   
    


 
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                            = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                        = "Vertical_Climb"   
    segment.analyses.extend(analyses.vertical_climb) 
    segment.altitude_start                             = 0.0  * Units.ft  
    segment.altitude_end                               = 40.  * Units.ft  
    segment.initial_battery_state_of_charge            = 1.0 
    segment.climb_rate                                 = 100. * Units['ft/min'] 

    # define flight dynamics to model  
    segment.flight_dynamics.force_z                        = True 

    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                            'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    
    mission.append_segment(segment)  
      
  
          
    # ------------------------------------------------------------------
    #  First Transition Segment
    # ------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                                           = "Vertical_Transition"  
    segment.analyses.extend( analyses.vertical_transition)   
    segment.air_speed_end                                 = 35 * Units['mph']     
    segment.acceleration                                  = 0.5

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True 
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                      = "low_speed_climb_transition" 
    segment.analyses.extend(analyses.climb_transition) 
    segment.climb_rate               = 500. * Units['ft/min'] 
    segment.air_speed_end            = 85.   * Units['mph'] 
    segment.altitude_start           = 40.0 * Units.ft    
    segment.altitude_end             = 100.0 * Units.ft

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True 
                                                                             
    mission.append_segment(segment)   
    
     
    # ------------------------------------------------------------------
    #  Second Transition Segment
    # ------------------------------------------------------------------ 
    segment                           = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                       = "high_speed_climb_transition"  
    segment.analyses.extend( analyses.climb_transition)   
    segment.air_speed_end             = 125.  * Units['mph']  
    segment.acceleration              = 9.81/5 

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                           = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                       = "Climb"  
    segment.analyses.extend(analyses.cruise) 
    segment.climb_rate                = 500. * Units['ft/min']
    segment.air_speed_start           = 125.   * Units['mph']
    segment.air_speed_end             = 130.  * Units['mph']  
    segment.altitude_start            = 100.0 * Units.ft   
    segment.altitude_end              = 2500.0 * Units.ft 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]  
        
    segment.assigned_control_variables.body_angle.active             = True
    
    mission.append_segment(segment)
    
    

    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                      = "Cruise"  
    segment.analyses.extend(analyses.cruise) 
    segment.altitude                 = 2500.0 * Units.ft
    segment.air_speed                = 130.  * Units['mph']   
    segment.distance                 = 30*Units.nmi
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #    Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                      = "Descent"  
    segment.analyses.extend(analyses.cruise)
    segment.climb_rate               = -300. * Units['ft/min']
    segment.air_speed_start          = 130.  * Units['mph'] 
    segment.air_speed_end            = 100.   * Units['mph'] 
    segment.altitude_start           = 2500.0 * Units.ft
    segment.altitude_end             = 100.0 * Units.ft

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
    
        
    mission.append_segment(segment)     

    # ------------------------------------------------------------------
    #   Reserve Climb Segment 
    # ------------------------------------------------------------------ 
    segment                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                      = "Reserve_Climb"   
    segment.analyses.extend(analyses.cruise) 
    segment.climb_rate               = 500. * Units['ft/min']
    segment.air_speed_start          = 100.   * Units['mph'] 
    segment.air_speed_end            = 120.  * Units['mph']  
    segment.altitude_start           = 100.0 * Units.ft 
    segment.altitude_end             = 1000.0 * Units.ft
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
             
    mission.append_segment(segment)      
 
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag                      = "Reserve_Cruise"  
    segment.analyses.extend(analyses.cruise)  
    segment.air_speed                = 120.  * Units['mph']  
    segment.distance                 = 10.*Units.nmi

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
    
        
    mission.append_segment(segment)     
 
    # ------------------------------------------------------------------
    #   Reserve Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Descent.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                      = "Reserve_Descent" 
    segment.analyses.extend(analyses.cruise)
    segment.descent_rate             = 300. * Units['ft/min']
    segment.air_speed_start          = 120.  * Units['mph'] 
    segment.air_speed_end            = 85.   * Units['mph']
    segment.altitude_start           = 1000.0 * Units.ft
    segment.altitude_end             = 100.0 * Units.ft

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True 
        
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Forth Transition Segment
    # ------------------------------------------------------------------ 
    segment                          = Segments.Descent.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                      = "Approach_Transition"   
    segment.analyses.extend(analyses.approach_transition)  
    segment.descent_rate             = 50.  * Units['ft/min'] 
    segment.air_speed_end            = 200. * Units['ft/min']   #20.   * Units['mph']   
    segment.altitude_end             = 40.0 * Units.ft

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]
    segment.assigned_control_variables.body_angle.active             = True
    
        
    mission.append_segment(segment)
    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Vertical Descent 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                                         = Segments.Vertical_Flight.Descent(base_segment)
    segment.tag                                                     = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_descent)               
    segment.altitude_start                                          = 100.0 * Units.ft   
    segment.altitude_end                                            = 0.   * Units.ft  
    segment.descent_rate                                            = 200. * Units['ft/min']  
                  
    # define flight dynamics to model              
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3','lift_rotor_propulsor_4',
                                                                             'lift_rotor_propulsor_5','lift_rotor_propulsor_6','lift_rotor_propulsor_7','lift_rotor_propulsor_8']]  
            
    mission.append_segment(segment)      
    
    return mission


def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  

# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------

def plot_results(results):
    # Plots fligh conditions 
    plot_flight_conditions(results) 
    
    # Plot arcraft trajectory
    plot_flight_trajectory(results)   

    plot_propulsor_throttles(results)
    
    # Plot Aircraft Electronics
    plot_battery_module_conditions(results) 
    plot_battery_temperature(results)
    plot_battery_cell_conditions(results) 
    plot_battery_module_C_rates(results)
    plot_battery_degradation(results) 
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
    plot_disc_and_power_loading(results)
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_propulsor_efficiencies(results)
    
      
    return  


if __name__ == '__main__': 
    main()    
    plt.show()
