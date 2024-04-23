# Vehicle.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     

import RCAIDE
from RCAIDE.Framework.Core import Units , Data
import numpy as np    
from RCAIDE.Plots.Performance.Mission_Plots                                   import *  
from RCAIDE.Plots.Geometry                                                    import * 
from RCAIDE.Methods.Propulsion                                                import design_propeller  

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup():

    '''
    This function defines the base vehicle including 
    1) center of gravity (either hard coded or use RCAIDE's built in function)
    2) mass moment of interita (optional)
    
    Key Notes:
    1) The wing that is intended to be the main must be given the tag "main wing". This wing will be used to append 
       a flap and an aileron 
    
    2) If present, the wing that is intended to be the horizontal stabilizer must be given the tag "horizontal_stabilizer" 
       This wing will be used to append an elevator 
    
    3) If present, The wing that is intended to be the  vertical stabilizer must be given the tag "vertical_stabilizer" 
       This wing will be used to append a rudder (optional)
    
    
    ''' 
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------ 
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Navion' 

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    vehicle.mass_properties.max_takeoff   = 1335.76 
    vehicle.mass_properties.takeoff       = 1335.76   
    vehicle.mass_properties.moments_of_inertia.tensor = np.array([[164627.7,0.0,0.0],[0.0,471262.4,0.0],[0.0,0.0,554518.7]])
    vehicle.mass_properties.center_of_gravity = [[2.239696797,0,-0.131189711 ]]
    vehicle.envelope.ultimate_load        = 5.7
    vehicle.envelope.limit_load           = 3.8
    vehicle.reference_area                = 17.112 
    vehicle.passengers                    = 2
      
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------    
    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.leading_edge              = 2.996 * Units.degrees 
    wing.thickness_to_chord               = 12
    wing.areas.reference                  = 17.112
    wing.chords.mean_aerodynamic          = 1.74 
    wing.taper                            = 0.54 
    wing.aspect_ratio                     = 6.04  
    wing.spans.projected                  = np.sqrt(wing.aspect_ratio*wing.areas.reference ) 
    wing.chords.root                      = wing.chords.mean_aerodynamic/( 2/3*(( 1 + wing.taper+ wing.taper**2 ) /( 1 + wing.taper )))  
    wing.chords.tip                       = wing.chords.root* wing.taper   
    wing.twists.root                      = 2 * Units.degrees  
    wing.twists.tip                       = -1 * Units.degrees   
    wing.dihedral                         = 7.5 * Units.degrees   
    wing.origin                           = [[1.652555594, 0.,-0.6006666]]
    wing.aerodynamic_center               = [1.852555594, 0., 6006666 ] # INCORRECT 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    
    tip_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = '../Airfoils/NACA_6410.txt' 
 
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    root_airfoil.coordinate_file          = '../Airfoils/NACA_4415.txt'  
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------       
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'horizontal_stabilizer'  
    wing.sweeps.leading_edge              = 6 * Units.degrees 
    wing.thickness_to_chord               = 12
    wing.areas.reference                  = 4  
    wing.taper                            = 0.67 
    wing.aspect_ratio                     = 4 
    wing.spans.projected                  = np.sqrt(wing.aspect_ratio*wing.areas.reference )
    wing.chords.root                      = 1.239465779 
    wing.chords.mean_aerodynamic          = wing.chords.root  *  2/3*(( 1 + wing.taper+ wing.taper**2 ) /( 1 + wing.taper ))
    wing.chords.tip                       = wing.chords.root* wing.taper   
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.origin                           = [[ 6.54518625 , 0., 0.203859697]]
    wing.aerodynamic_center               = [ 6.545186254 + 0.25*wing.spans.projected, 0., 0.203859697]
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
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'vertical_stabilizer'   
    wing.sweeps.leading_edge              = 20 * Units.degrees 
    wing.thickness_to_chord               = 12.5
    wing.areas.reference                  = 1.163   
    wing.spans.projected                  = 1.481640283 
    wing.chords.root                      = 1.217672543
    wing.chords.tip                       = 0.587043035
    wing.aspect_ratio                     = wing.spans.projected**2/ wing.areas.reference 
    wing.taper                            = wing.chords.tip/wing.chords.root
    wing.chords.mean_aerodynamic          = wing.chords.root  *  2/3*(( 1 + wing.taper+ wing.taper**2 ) /( 1 + wing.taper ))
    wing.chords.tip                       = wing.chords.root 
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.origin                           = [[ 7.127369987, 0., 0.303750948]]
    wing.aerodynamic_center               = [ 7.49778005775, 0., 0.67416101875] 
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
    fuselage = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.seats_abreast                      = 2
    fuselage.lengths.total                      = 8.349950916 
    fuselage.width                              = 1.22028016 
    fuselage.heights.maximum                    = 1.634415138  
    fuselage.areas.wetted                       = 12. # ESTIMATED 
    fuselage.areas.front_projected              = fuselage.width*fuselage.heights.maximum
    fuselage.effective_diameter                 = 1.22028016 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.529255748
    segment.width                               = 0.575603849
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  =  0.028527593
    segment.percent_z_location                  =  0
    segment.height                              =  0.737072721
    segment.width                               =  0.921265952 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.187342754 
    segment.percent_z_location                  = 0 
    segment.height                              = 1.174231852 
    segment.width                               = 1.196956212
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.242034847 
    segment.percent_z_location                  = 0.011503528 
    segment.height                              = 1.450221906 
    segment.width                               = 1.173932059 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.296715183 
    segment.percent_z_location                  = 0.015984303 
    segment.height                              = 1.634415138 
    segment.width                               = 1.22028016 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.510275342 
    segment.percent_z_location                  = -0.005
    segment.height                              = 1.082135236 
    segment.width                               = 1.013062774 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.833284347 
    segment.percent_z_location                  = 0.014138855 
    segment.height                              = 0.621652157 
    segment.width                               = 0.414134978
    fuselage.Segments.append(segment)
 
    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.018978667 
    segment.height                              = 0.092096616 
    segment.width                               = 0.046048308 
    fuselage.Segments.append(segment)
    
    # add to vehicle
    vehicle.append_component(fuselage)
 
    # ------------------------------------------------------------------
    #   Fuel
    # ------------------------------------------------------------------    
    # define fuel weight needed to size fuel system
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline()
    fuel.mass_properties                        = RCAIDE.Library.Components.Mass_Properties() 
    fuel.number_of_tanks                        = 1.
    fuel.origin                                 = wing.origin
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density #all of the fuel volume is internal
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.mass_properties.mass                   = 319 *Units.lbs
    vehicle.fuel                                = fuel

    # ------------------------------------------------------------------
    #   Piston Propeller Network
    # ------------------------------------------------------------------    
    
    # build network
    net                                         = RCAIDE.Framework.Networks.Internal_Combustion_Propeller()
    net.tag                                     = 'internal_combustion'
    net.number_of_engines                       = 1.
    net.identical_propellers                    = True
                                                
    # the engine                    
    engine                                  = RCAIDE.Library.Components.Propulsors.Converters.Internal_Combustion_Engine()
    engine.sea_level_power                  = 180. * Units.horsepower
    engine.flat_rate_altitude               = 0.0
    engine.rated_speed                      = 2700. * Units.rpm
    engine.power_specific_fuel_consumption  = 0.52 
    net.engines.append(engine)
    
    # the prop
    prop = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
    prop.number_of_blades                      = 2.0
    prop.tip_radius                            = 2.14/2
    prop.hub_radius                            = 8.     * Units.inches
    prop.cruise.design_Cl                      = 0.8
    prop.cruise.design_altitude                = 12000. * Units.feet
    prop.cruise.design_power                   = .64 * 180. * Units.horsepower
    prop.cruise.design_freestream_velocity     = 119.   * Units.knots
    prop.cruise.design_angular_velocity        = 2650.  * Units.rpm
    prop.variable_pitch                        = True
    prop.origin                                = [[-0.1, 0, 0]] 
    prop.airfoil_geometry                      = '../Airfoils/NACA_4412.txt'
    prop.airfoil_polars                        = ['../Airfoils/Polars/NACA_4412_polar_Re_50000.txt' ,
                                                   '../Airfoils/Polars/NACA_4412_polar_Re_100000.txt' ,
                                                   '../Airfoils/Polars/NACA_4412_polar_Re_200000.txt' ,
                                                   '../Airfoils/Polars/NACA_4412_polar_Re_500000.txt' ,
                                                   '../Airfoils/Polars/NACA_4412_polar_Re_1000000.txt' ]
               
    prop.airfoil_polar_stations                = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]      
    prop                                       = design_propeller(prop)   
    
    net.rotors.append(prop)
     
    
    # add the network to the vehicle
    vehicle.append_component(net) 

    #find uninstalled avionics weight
    Wuav                                        = 2. * Units.lbs
    avionics                                    = RCAIDE.Library.Components.Systems.Avionics()
    avionics.mass_properties.uninstalled        = Wuav
    vehicle.avionics                            = avionics  
    
    
    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------  
    #plot_3d_vehicle(vehicle) 
    #plt.show()  
    

    return vehicle 

def stick_fixed_stability_setup(): 
    vehicle  = vehicle_setup()   
    configs  = stick_fixed_stability_configs_setup(vehicle) 
    return configs 
 
 
def elevator_sizing_setup(vehicle):   
    hs_wing                        = vehicle.wings.horizontal_stabilizer 
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.1
    elevator.span_fraction_end     = 0.9
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    hs_wing.append_control_surface(elevator)     
    configs                        = elevator_sizing_configs_setup(vehicle) 
    return configs

def aileron_rudder_sizing_setup(vehicle):    
    mw_wing                       = vehicle.wings.main_wing 
    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.9 
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.2
    mw_wing.append_control_surface(aileron) 
    
    if vehicle.rudder_flag:
        vs_wing                      = vehicle.wings.vertical_stabilizer 
        rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
        rudder.tag                   = 'rudder'
        rudder.span_fraction_start   = 0.2
        rudder.span_fraction_end     = 0.8
        rudder.deflection            = 0.0  * Units.deg
        rudder.chord_fraction        = 0.2
        vs_wing.append_control_surface(rudder) 
    
    configs  = aileron_rudder_sizing_configs_setup(vehicle) 
    return configs 
 
def flap_sizing_setup(vehicle):   
    mw_wing                       = vehicle.wings.main_wing
    flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.5
    flap.deflection               = 0.0 * Units.degrees 
    flap.chord_fraction           = 0.20
    mw_wing.append_control_surface(flap)   
    configs                       = flap_sizing_configs_setup(vehicle) 
    return configs 

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def stick_fixed_stability_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle) 
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'stick_fixed_cruise'
    configs.append(config) 
    return configs  

def elevator_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle) 
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'elevator_sizing'   
    configs.append(config)     
    return configs

def aileron_rudder_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)  
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'aileron_rudder_sizing'   
    configs.append(config)   
    return configs  

def flap_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container()  
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)  
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'flap_sizing'   
    configs.append(config)   
    return configs
