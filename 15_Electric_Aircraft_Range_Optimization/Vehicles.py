# Vehicle.py 
 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy as np    
import os
import pickle
from copy import deepcopy 

import RCAIDE 
from RCAIDE.Core                                                              import Data , Units 
from RCAIDE.Energy.Networks.Battery_Electric_Rotor                            import Battery_Electric_Rotor
from RCAIDE.Methods.Propulsion                                                import design_propeller 
from RCAIDE.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration  
from RCAIDE.Methods.Propulsion.electric_motor_sizing                          import size_optimal_motor
from RCAIDE.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from RCAIDE.Methods.Weights.Buildups.eVTOL.empty                              import empty  
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_segmented_planform import wing_segmented_planform
from RCAIDE.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
from RCAIDE.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity
 
# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------    

def setup():
    
    base_vehicle = base_vehicle_setup()
    configs      = configs_setup(base_vehicle)
    
    return configs

def base_vehicle_setup() :
     
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------     
    vehicle                               = RCAIDE.Vehicle() 
    vehicle.tag                           = 'electric_CTOL'

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    vehicle.mass_properties.max_takeoff   = 1900  # kg 
    vehicle.mass_properties.takeoff       = 1900  # kg 
    vehicle.mass_properties.max_zero_fuel = 1900  # kg 
    vehicle.envelope.ultimate_load        = 5.7
    vehicle.envelope.limit_load           = 3.8 
    vehicle.reference_area                = 14.76
    vehicle.passengers                    = 6
    vehicle.systems.control               = "fully powered"
    vehicle.systems.accessories           = "commuter"    
    
    cruise_speed                          = 175.*Units['mph']    
    altitude                              = 2500. * Units.ft
    atmo                                  = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    freestream                            = atmo.compute_values (0.)
    freestream0                           = atmo.compute_values (altitude)
    mach_number                           = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure       = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.design_mach_number            =  mach_number
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------    
    wing                                  = RCAIDE.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 14.76
    wing.spans.projected                  = 11.2 
    wing.chords.root                      = 1.2 #1.46
    wing.chords.tip                       = 0.9 # 0.92
    wing.chords.mean_aerodynamic          = 1.19
    wing.taper                            = wing.chords.root/wing.chords.tip 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 4.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[2.93, 0., 1.01]]
    wing.aerodynamic_center               = [3., 0., 1.01] 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    airfoil                               = RCAIDE.Components.Airfoils.Airfoil()
    airfoil.coordinate_file               = '../Airfoils/NACA_63_412.txt'
    
    cg_x = wing.origin[0][0] + 0.25*wing.chords.mean_aerodynamic
    cg_z = wing.origin[0][2] - 0.2*wing.chords.mean_aerodynamic
    vehicle.mass_properties.center_of_gravity = [[cg_x,   0.  ,  cg_z ]]  # SOURCE: Design and aerodynamic analysis of a twin-engine commuter aircraft

    # Wing Segments
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'inboard'
    segment.percent_span_location         = 0.0 
    segment.twist                         = 4. * Units.degrees   
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0.  
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)

    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'outboard'
    segment.percent_span_location         = 0.5438
    segment.twist                         = 3.0* Units.degrees 
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0. 
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)
    
    # Wing Segments
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'winglet'
    segment.percent_span_location         = 0.98
    segment.twist                         = 0.  * Units.degrees 
    segment.root_chord_percent            = 0.60
    segment.dihedral_outboard             = 75. * Units.degrees 
    segment.sweeps.quarter_chord          = 82. * Units.degrees 
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(airfoil)
    wing.append_segment(segment) 

    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees 
    segment.root_chord_percent            = 0.12
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)    
    
    # Fill out more segment properties automatically 
    wing_segmented_planform(wing, overwrite_reference = True ) 
    wing = segment_properties(wing)
    vehicle.reference_area                = wing.areas.reference   
    
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------       
    wing                                  = RCAIDE.Components.Wings.Wing()
    wing.tag                              = 'horizontal_stabilizer' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 2.540 
    wing.spans.projected                  = 3.3  * Units.meter 
    wing.sweeps.quarter_chord             = 0 * Units.deg 
    wing.chords.root                      = 0.769 * Units.meter 
    wing.chords.tip                       = 0.769 * Units.meter 
    wing.chords.mean_aerodynamic          = 0.769 * Units.meter  
    wing.taper                            = 1. 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[7.7, 0., 0.25]]
    wing.aerodynamic_center               = [7.8, 0., 0.25] 
    wing.vertical                         = False
    wing.winglet_fraction                 = 0.0  
    wing.symmetric                        = True
    wing.high_lift                        = False 
    wing.dynamic_pressure_ratio           = 0.9

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                                  = RCAIDE.Components.Wings.Wing()
    wing.tag                              = 'vertical_stabilizer'     
    wing.sweeps.quarter_chord             = 25. * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 2.258 * Units['meters**2']  
    wing.spans.projected                  = 1.854   * Units.meter  
    wing.chords.root                      = 1.6764 * Units.meter 
    wing.chords.tip                       = 0.6858 * Units.meter 
    wing.chords.mean_aerodynamic          = 1.21   * Units.meter 
    wing.taper                            = wing.chords.tip/wing.chords.root 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[6.75 ,0,  0.514]]
    wing.aerodynamic_center               = [0.508 ,0,0]  
    wing.vertical                         = True 
    wing.symmetric                        = False
    wing.t_tail                           = False
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    fuselage = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.seats_abreast                      = 2.
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2.
    fuselage.lengths.nose                       = 60.  * Units.inches
    fuselage.lengths.tail                       = 161. * Units.inches
    fuselage.lengths.cabin                      = 105. * Units.inches
    fuselage.lengths.total                      = 332.2* Units.inches
    fuselage.lengths.fore_space                 = 0.
    fuselage.lengths.aft_space                  = 0.
    fuselage.width                              = 42. * Units.inches
    fuselage.heights.maximum                    = 62. * Units.inches
    fuselage.heights.at_quarter_length          = 62. * Units.inches
    fuselage.heights.at_three_quarters_length   = 62. * Units.inches
    fuselage.heights.at_wing_root_quarter_chord = 23. * Units.inches
    fuselage.areas.side_projected               = 8000.  * Units.inches**2.
    fuselage.areas.wetted                       = 30000. * Units.inches**2.
    fuselage.areas.front_projected              = 42.* 62. * Units.inches**2.
    fuselage.effective_diameter                 = 50. * Units.inches 

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.0
    segment.width                               = 0.0
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.007279116466
    segment.percent_z_location                  = 0.001253114988
    segment.height                              = 0.1669064748
    segment.width                               = 0.2780205877
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.01941097724
    segment.percent_z_location                  = 0.001216095397
    segment.height                              = 0.3129496403
    segment.width                               = 0.4365777215
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.06308567604
    segment.percent_z_location                  = 0.007395489231
    segment.height                              = 0.5841726619
    segment.width                               = 0.6735119903
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.1653761217
    segment.percent_z_location                  = 0.02891281352
    segment.height                              = 1.064028777
    segment.width                               = 1.067200529
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.2426372155
    segment.percent_z_location                  = 0.04214148761
    segment.height                              = 1.293766653
    segment.width                               = 1.183058255
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.2960174029
    segment.percent_z_location                  = 0.04705241831
    segment.height                              = 1.377026712
    segment.width                               = 1.181540054
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.3809404284
    segment.percent_z_location                  = 0.05313580461
    segment.height                              = 1.439568345
    segment.width                               = 1.178218989
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.5046854083
    segment.percent_z_location                  = 0.04655492473
    segment.height                              = 1.29352518
    segment.width                               = 1.054390707
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  = 0.6454149933
    segment.percent_z_location                  = 0.03741966266
    segment.height                              = 0.8971223022
    segment.width                               = 0.8501926505
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.985107095
    segment.percent_z_location                  = 0.04540283436
    segment.height                              = 0.2920863309
    segment.width                               = 0.2012565415
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.04787575562
    segment.height                              = 0.1251798561
    segment.width                               = 0.1206021048
    fuselage.Segments.append(segment)

    # add to vehicle
    vehicle.append_component(fuselage)

    # ------------------------------------------------------------------
    #   Nacelles
    # ------------------------------------------------------------------ 
    nacelle                = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.tag            = 'nacelle_1'
    nacelle.length         = 2
    nacelle.diameter       = 42 * Units.inches
    nacelle.areas.wetted   = 0.01*(2*np.pi*0.01/2)
    nacelle.origin         = [[2.5,2.5,1.0]]
    nacelle.flow_through   = False  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.1  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.3  
    nac_segment.height             = 0.52
    nac_segment.width              = 0.7
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.5  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.7 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.6
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.9 
    nac_segment.height             = 0.3
    nac_segment.width              = 0.5
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 1.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)    
    
    vehicle.append_component(nacelle)  

    nacelle_2          = deepcopy(nacelle)
    nacelle_2.tag      = 'nacelle_2'
    nacelle_2.origin   = [[2.5,-2.5,1.0]]
    vehicle.append_component(nacelle_2)    
    
    #---------------------------------------------------------------------------------------------
    # BUILD ENERGY NETWORK
    #--------------------------------------------------------------------------------------------- 
    net                     = Battery_Electric_Rotor() 
    net.rotor_group_indexes = [0,0]
    net.motor_group_indexes = [0,0]  
    net.esc_group_indexes   = [0,0]    

    # Component 1: ESC
    esc_1                     = RCAIDE.Energy.Propulsion.Modulators.Electronic_Speed_Controller()
    esc_1.efficiency          = 0.95 
    esc_1.tag                 = 'esc_1'
    net.electronic_speed_controllers.append(esc_1)     

    esc_2                     = RCAIDE.Energy.Propulsion.Modulators.Electronic_Speed_Controller()
    esc_2.efficiency          = 0.95 
    esc_2.tag                 = 'esc_2'
    net.electronic_speed_controllers.append(esc_2)      
     
    
    # Component 2: Propeller 
    propeller                                   = RCAIDE.Energy.Propulsion.Converters.Propeller()
    propeller.number_of_blades                  = 3 
    propeller.tag                               = 'propeller_1'  
    propeller.tip_radius                        = 1.72/2   
    propeller.hub_radius                        = 10.     * Units.inches 
    propeller.cruise.design_freestream_velocity = 175.*Units['mph']   
    propeller.cruise.design_angular_velocity    = 2700. * Units.rpm 
    propeller.cruise.design_Cl                  = 0.7 
    propeller.cruise.design_altitude            = 2500. * Units.feet 
    propeller.cruise.design_thrust              = 2000   
    propeller.rotation                          = -1 
    propeller.variable_pitch                    = True  
    propeller.origin                            = [[2.,2.5,0.95]]
    airfoil                                     = RCAIDE.Components.Airfoils.Airfoil() 
    ospath                                      = os.path.abspath(__file__)
    separator                                   = os.path.sep
    rel_path                                    = os.path.dirname(ospath) + separator  
    airfoil.coordinate_file                     =  rel_path + '../Airfoils/NACA_4412.txt'
    airfoil.polar_files                         = [rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_50000.txt' ,
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_100000.txt' ,
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_200000.txt' ,
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_500000.txt' ,
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_1000000.txt',
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_3500000.txt',
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_5000000.txt',
                                                   rel_path + '../Airfoils/Polars/NACA_4412_polar_Re_7500000.txt' ]
    propeller.append_airfoil(airfoil)          
    propeller.airfoil_polar_stations            = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    propeller                                   = design_propeller(propeller) 
    net.rotors.append(propeller)   
        
    # append propeller on left hans side of wing
    propeller_left                               = deepcopy(propeller)
    propeller_left.tag                           = 'propeller_2' 
    propeller_left.origin                        = [[2.,-2.5,0.95]]
    propeller_left.rotation                      = 1 
    net.rotors.append(propeller_left) 
                 
    # Component 3: Battery       
    bat                                                    = RCAIDE.Energy.Sources.Batteries.Lithium_Ion_NMC() 
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 100
    initialize_from_circuit_configuration(bat)  
    bat.module_config.number_of_modules                    = 14  
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module_config.voltage                              = bat.pack.max_voltage/bat.module_config.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat.module.geometrtic_configuration.normal_count       = 24
    bat.module.geometrtic_configuration.parallel_count     = 40
    net.battery                                            = bat    
    net.voltage                                            = bat.pack.max_voltage
    
    # Component 4: Miscellaneous Systems
    sys                       = RCAIDE.Components.Systems.System()
    sys.mass_properties.mass  = 5 # kg
 
    # Component 5: Motor   
    propeller_motor                         = RCAIDE.Energy.Propulsion.Converters.Motor()
    propeller_motor.efficiency              = 0.95
    propeller_motor.gearbox_efficiency      = 1.
    propeller_motor.origin                  = [[2.,  2.5, 0.784]]
    propeller_motor.nominal_voltage         = bat.pack.max_voltage*0.8 
    propeller_motor.no_load_current         = 0.1 
    propeller_motor.rotor_radius            = propeller.tip_radius
    propeller_motor.design_torque           = propeller.cruise.design_torque
    propeller_motor.angular_velocity        = propeller.cruise.design_angular_velocity/propeller_motor.gear_ratio
    propeller_motor                         = size_optimal_motor(propeller_motor)  
    propeller_motor.mass_properties.mass    = nasa_motor(propeller_motor.design_torque)   
    net.motors.append(propeller_motor)   
    
    # append left motor 
    propeller_motor_left                        = deepcopy(propeller_motor)
    propeller_motor_left.origin                 = [[2., -2.5, 0.784]] 
    net.motors.append(propeller_motor_left) 

    # Component 6: Payload
    payload                                     = RCAIDE.Energy.Peripherals.Payload()
    payload.power_draw                          = 10. # Watts
    payload.mass_properties.mass                = 1.0 * Units.kg
    net.payload                                 = payload
               
    # Component 7: Avionics               
    avionics                                    = RCAIDE.Energy.Peripherals.Avionics()
    avionics.power_draw                         = 20. # Watts
    net.avionics                                = avionics

    motor_origins                                         = np.array([[2., -2.5, 0.784],[2., 2.5, 0.784]])
    vehicle.wings['main_wing'].motor_spanwise_locations   = motor_origins[:,1]/ vehicle.wings['main_wing'].spans.projected
    
    # add the solar network to the vehicle
    vehicle.append_component(net)
     
    settings = Data()
    converge_evtol_weight(vehicle,settings,contingency_factor = 1.1) 
    breakdown = empty(vehicle,settings,contingency_factor     = 1.1 ) 
    
    vehicle.weight_breakdown  = breakdown
    compute_component_centers_of_gravity(vehicle)
    vehicle.center_of_gravity()  
    
    return vehicle 
     
 

# ---------------------------------------------------------------------
#   Define the Configurations
# --------------------------------------------------------------------- 
def configs_setup(vehicle): 
    
    configs = RCAIDE.Components.Configs.Config.Container()

    base_config = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)    
     
    return configs 

 