# tut_mission_B737_AVL.py
# 
# Created:  Mar 2018, RCAIDE Team

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE Imports
import RCAIDE 
from RCAIDE.Core import Units
from RCAIDE.Visualization  import *    
from RCAIDE.Methods.Propulsion.turbofan_sizing import turbofan_sizing 
import numpy as np  
import matplotlib.pyplot as plt
from copy import deepcopy

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
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Analyses.Vehicle()
 
    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.AVL()
    aerodynamics.process.compute.lift.inviscid.settings.filenames.avl_bin_name = 'CHANGE ME TO YOUR DIRECTORY' # eg. '/Users/matthewclarke/Documents/AVL/avl3.35'
    #aerodynamics.settings.number_spanwise_vortices  = 5
    #aerodynamics.settings.number_chordwise_vortices = 3
    aerodynamics.geometry = vehicle
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = RCAIDE.Analyses.Stability.AVL()
    stability.settings.filenames.avl_bin_name = 'CHANGE ME TO YOUR DIRECTORY' # eg. '/Users/matthewclarke/Documents/AVL/avl3.35'
    
    stability.geometry = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.networks = vehicle.networks
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses    

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():
    """This is the full physical definition of the vehicle, and is designed to be independent of the
    analyses that are selected."""
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Boeing_737-800'    
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------     
    vehicle.mass_properties.max_takeoff               = 79015.8 * Units.kilogram  
    vehicle.mass_properties.takeoff                   = 79015.8 * Units.kilogram    
    vehicle.mass_properties.operating_empty           = 62746.4 * Units.kilogram  
    vehicle.mass_properties.max_zero_fuel             = 62732.0 * Units.kilogram 
    vehicle.mass_properties.cargo                     = 10000.  * Units.kilogram    
    vehicle.envelope.ultimate_load = 3.75
    vehicle.envelope.limit_load    = 2.5

    # Vehicle level parameters 
    vehicle.reference_area         = 124.862 * Units['meters**2']    
    vehicle.passengers             = 170
    vehicle.systems.control        = "fully powered" 
    vehicle.systems.accessories    = "medium range"

    # ------------------------------------------------------------------        
    #  Landing Gear
    # ------------------------------------------------------------------  
    landing_gear = RCAIDE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag = "main_landing_gear"
    
    landing_gear.main_tire_diameter = 1.12000 * Units.m
    landing_gear.nose_tire_diameter = 0.6858 * Units.m
    landing_gear.main_strut_length  = 1.8 * Units.m
    landing_gear.nose_strut_length  = 1.3 * Units.m
    landing_gear.main_units  = 2    # Number of main landing gear
    landing_gear.nose_units  = 1    # Number of nose landing gear
    landing_gear.main_wheels = 2    # Number of wheels on the main landing gear
    landing_gear.nose_wheels = 2    # Number of wheels on the nose landing gear      
    vehicle.landing_gear = landing_gear

    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------         
    
    wing = RCAIDE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.aspect_ratio            = 10.18 
    wing.sweeps.quarter_chord    = 25 * Units.deg
    wing.thickness_to_chord      = 0.1
    wing.taper                   = 0.1
    wing.spans.projected         = 34.32 * Units.meter
    wing.chords.root             = 7.760 * Units.meter
    wing.chords.tip              = 0.782 * Units.meter
    wing.chords.mean_aerodynamic = 4.235 * Units.meter
    wing.areas.reference         = 124.862 * Units['meters**2']  
    wing.twists.root             = 4.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees
    wing.origin                  = [[13.61, 0, -1.27]] * Units.meter
    wing.vertical                = False
    wing.symmetric               = True 
    wing.high_lift               = True 
    wing.dynamic_pressure_ratio  = 1.0
    
    # ------------------------------------------------------------------
    #   Main Wing Control Surfaces
    # ------------------------------------------------------------------ 
    flap                       = RCAIDE.Components.Wings.Control_Surfaces.Flap() 
    flap.tag                   = 'flap' 
    flap.span_fraction_start   = 0.20 
    flap.span_fraction_end     = 0.70   
    flap.deflection            = 0.0 * Units.degrees 
    flap.configuration_type    = 'double_slotted'
    flap.chord_fraction        = 0.30   
    wing.append_control_surface(flap)   
        
    slat                       = RCAIDE.Components.Wings.Control_Surfaces.Slat() 
    slat.tag                   = 'slat' 
    slat.span_fraction_start   = 0.324 
    slat.span_fraction_end     = 0.963     
    slat.deflection            = 0.0 * Units.degrees
    slat.chord_fraction        = 0.1  	 
    wing.append_control_surface(slat)  
        
    aileron                       = RCAIDE.Components.Wings.Control_Surfaces.Aileron() 
    aileron.tag                   = 'aileron' 
    aileron.span_fraction_start   = 0.7 
    aileron.span_fraction_end     = 0.963 
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16    
    wing.append_control_surface(aileron)    
    
    # Add to vehicle
    vehicle.append_component(wing)    

    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------        
    
    wing = RCAIDE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'
    
    wing.aspect_ratio            = 6.16     
    wing.sweeps.quarter_chord    = 40.0 * Units.deg
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.2
    wing.spans.projected         = 14.2 * Units.meter
    wing.chords.root             = 4.7  * Units.meter
    wing.chords.tip              = 0.955 * Units.meter
    wing.chords.mean_aerodynamic = 3.0  * Units.meter
    wing.areas.reference         = 32.488   * Units['meters**2']  
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees  
    wing.origin                  = [[32.83 * Units.meter, 0 , 1.14 * Units.meter]]
    wing.vertical                = False 
    wing.symmetric               = True
    wing.dynamic_pressure_ratio  = 0.9  
    
    # Add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'    

    wing.aspect_ratio            = 1.91
    wing.sweeps.quarter_chord    = 25. * Units.deg
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.25
    wing.spans.projected         = 7.777 * Units.meter
    wing.chords.root             = 8.19  * Units.meter
    wing.chords.tip              = 0.95  * Units.meter
    wing.chords.mean_aerodynamic = 4.0   * Units.meter
    wing.areas.reference         = 27.316 * Units['meters**2']  
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [[28.79 * Units.meter, 0, 1.54 * Units.meter]] # meters
    wing.vertical                = True 
    wing.symmetric               = False 
    wing.t_tail                  = False
    wing.dynamic_pressure_ratio  = 1.0
        
    # Add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage' 
    fuselage.number_coach_seats    = vehicle.passengers 
    fuselage.seats_abreast         = 6
    fuselage.seat_pitch            = 1     * Units.meter 
    fuselage.fineness.nose         = 1.6
    fuselage.fineness.tail         = 2. 
    fuselage.lengths.nose          = 6.4   * Units.meter
    fuselage.lengths.tail          = 8.0   * Units.meter
    fuselage.lengths.total         = 38.02 * Units.meter 
    fuselage.lengths.fore_space    = 6.    * Units.meter
    fuselage.lengths.aft_space     = 5.    * Units.meter
    fuselage.width                 = 3.74  * Units.meter
    fuselage.heights.maximum       = 3.74  * Units.meter
    fuselage.effective_diameter    = 3.74     * Units.meter
    fuselage.areas.side_projected  = 142.1948 * Units['meters**2'] 
    fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
    fuselage.areas.front_projected = 12.57    * Units['meters**2']  
    fuselage.differential_pressure = 5.0e4 * Units.pascal 
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter
    
    # add to vehicle
    vehicle.append_component(fuselage)
    
    
    # ------------------------------------------------------------------
    #   Nacelles
    # ------------------------------------------------------------------ 
    nacelle                       = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.tag                   = 'nacelle_1'
    nacelle.length                = 2.71
    nacelle.inlet_diameter        = 1.90
    nacelle.diameter              = 2.05
    nacelle.areas.wetted          = 1.1*np.pi*nacelle.diameter*nacelle.length
    nacelle.origin                = [[13.72, -4.86,-1.9]]
    nacelle.flow_through          = True  
    nacelle_airfoil               = RCAIDE.Components.Airfoils.Airfoil() 
    nacelle_airfoil.naca_4_series_airfoil = '2410'
    nacelle.append_airfoil(nacelle_airfoil)

    nacelle_2                     = deepcopy(nacelle)
    nacelle_2.tag                 = 'nacelle_2'
    nacelle_2.origin              = [[13.72, 4.86,-1.9]]
    
    vehicle.append_component(nacelle)  
    vehicle.append_component(nacelle_2)     
         
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Turbofan Network
    #------------------------------------------------------------------------------------------------------------------------------------   
    #initialize the gas turbine network
    net                   = RCAIDE.Energy.Networks.Turbofan_Engine() 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Energy.Distributors.Fuel_Line() 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Fuel
    #------------------------------------------------------------------------------------------------------------------------------------  
    # fuel tank
    fuel_tank                                   = RCAIDE.Energy.Storages.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin 
    
    # fuel 
    fuel                                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline()   
    fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel  
    fuel_line.fuel_tanks.append(fuel_tank) 
     
    turbofan                                    = RCAIDE.Energy.Converters.Turbofan() 
    turbofan.tag                                = 'turbofan'
    turbofan.origin                             = [[12.0,4.38,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 5.4   
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 37278.0* Units.N/2 
             
    # fan                
    fan                                         = RCAIDE.Energy.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        
                   
    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Attributes.Gases.Air() 
    ram                                         = RCAIDE.Energy.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 
          
    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle


    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Energy.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 1.9   
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Energy.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.0    
    turbofan.high_pressure_compressor              = high_pressure_compressor

    # low pressure turbine  
    low_pressure_turbine                        = RCAIDE.Energy.Converters.Turbine()   
    low_pressure_turbine.tag                    ='lpt'
    low_pressure_turbine.mechanical_efficiency  = 0.99
    low_pressure_turbine.polytropic_efficiency  = 0.93 
    turbofan.low_pressure_turbine               = low_pressure_turbine

    # high pressure turbine  
    high_pressure_turbine                       = RCAIDE.Energy.Converters.Turbine()   
    high_pressure_turbine.tag                   ='hpt'
    high_pressure_turbine.mechanical_efficiency = 0.99
    high_pressure_turbine.polytropic_efficiency = 0.93 
    turbofan.high_pressure_turbine              = high_pressure_turbine 

    # combustor  
    combustor                           = RCAIDE.Energy.Converters.Combustor()   
    combustor.tag                       = 'Comb'
    combustor.efficiency                = 0.99 
    combustor.alphac                    = 1.0     
    combustor.turbine_inlet_temperature = 1500
    combustor.pressure_ratio            = 0.95
    combustor.fuel_data                 = RCAIDE.Attributes.Propellants.Jet_A()  
    turbofan.combustor                  = combustor

    # core nozzle
    core_nozzle                       = RCAIDE.Energy.Converters.Expansion_Nozzle()   
    core_nozzle.tag                   = 'core nozzle'
    core_nozzle.polytropic_efficiency = 0.95
    core_nozzle.pressure_ratio        = 0.99  
    turbofan.core_nozzle              = core_nozzle

    # fan nozzle
    fan_nozzle                       = RCAIDE.Energy.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                   = 'fan nozzle'
    fan_nozzle.polytropic_efficiency = 0.95
    fan_nozzle.pressure_ratio        = 0.99 
    turbofan.fan_nozzle              = fan_nozzle 
    
    #design turbofan
    design_turbofan(turbofan) 
    
    fuel_line.turbofans.append(turbofan)

    turbofan_2                      = deepcopy(turbofan)
    turbofan_2.tag                  = 'propeller_2' 
    turbofan_2.origin               = [[12.0,-4.38,-1.1]] 
    fuel_line.turbofans.append(turbofan_2)    

    net.fuel_lines.append(fuel_line)    
    vehicle.append_energy_network(net)     
        
    return vehicle


# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = RCAIDE.Components.Configs.Config.Container() 
    base_config = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)

    return configs 

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):

   
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission' 
    Segments = RCAIDE.Analyses.Mission.Segments  
    base_segment = Segments.Segment()


    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 3.0   * Units.km
    segment.air_speed      = 125.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 190.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end = 10.5   * Units.km
    segment.air_speed    = 226.0  * Units['m/s']
    segment.climb_rate   = 3.0    * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude  = 10.668 * Units.km  
    segment.air_speed = 230.412 * Units['m/s']
    segment.distance  = 1000 * Units.nmi 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_start = 10.5 * Units.km 
    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 220.0 * Units['m/s']
    segment.descent_rate   = 4.5   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2" 
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end = 6.0   * Units.km
    segment.air_speed    = 195.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"  
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end = 4.0   * Units.km
    segment.air_speed    = 170.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Fourth Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_4" 
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end = 2.0   * Units.km
    segment.air_speed    = 150.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)



    # ------------------------------------------------------------------
    #   Fifth Descent Segment:Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_5" 
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 145.0 * Units['m/s']
    segment.descent_rate = 3.0   * Units['m/s'] 
    segment = analyses.base.energy.networks.turbofan_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(mission):


    missions     = RCAIDE.Analyses.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions  

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results):

      # Plot Flight Conditions 
    plot_flight_conditions(results)
    
    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results)
    
    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results)
    
    # Plot Static Stability Coefficients 
    plot_stability_coefficients(results)    
    
    # Drag Components
    plot_drag_components(results)
    
    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results)
    
    # Plot Velocities 
    plot_aircraft_velocities(results)  
        
    return

if __name__ == '__main__': 
    main()
    plt.show()