# tut_Concorde.py
 

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from RCAIDE.Core import Units   
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform      import segment_properties   
from RCAIDE.Methods.Propulsion                             import design_turbofan
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform      import wing_planform, segment_properties
from RCAIDE.Methods.Propulsion                             import design_turbojet
from RCAIDE.Visualization                                  import *     

# python imports 
import numpy as np  
from copy import deepcopy
import matplotlib.pyplot as plt  
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
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Supersonic_Zero()
    aerodynamics.geometry = vehicle    
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    aerodynamics.settings.span_efficiency            = .8
    
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    
    # Not yet available for this configuration

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.networks = vehicle.networks #what is called throughout the mission (at every time step))
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

    # done!
    return analyses    

def vehicle_setup():

    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Concorde'    
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    

    # mass properties
    vehicle.mass_properties.max_takeoff     = 185000. * Units.kilogram   
    vehicle.mass_properties.operating_empty = 78700.  * Units.kilogram   
    vehicle.mass_properties.takeoff         = 185000. * Units.kilogram   
    vehicle.mass_properties.cargo           = 1000.   * Units.kilogram   
        
    # envelope properties
    vehicle.envelope.ultimate_load = 3.75
    vehicle.envelope.limit_load    = 2.5

    # basic parameters
    vehicle.reference_area               = 358.25      
    vehicle.passengers                   = 100
    vehicle.systems.control              = "fully powered" 
    vehicle.systems.accessories          = "long range"
    vehicle.maximum_cross_sectional_area = 13.9
    vehicle.total_length                 = 61.66
    
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------        
    
    wing = RCAIDE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.aspect_ratio            = 1.83
    wing.sweeps.quarter_chord    = 59.5 * Units.deg
    wing.thickness_to_chord      = 0.03
    wing.taper                   = 0.
    wing.spans.projected         = 25.6 * Units.meter
    wing.chords.root             = 33.8 * Units.meter
    wing.total_length            = 33.8 * Units.meter
    wing.chords.tip              = 1.1  * Units.meter
    wing.chords.mean_aerodynamic = 18.4 * Units.meter
    wing.areas.reference         = 358.25 * Units['meter**2']  
    wing.areas.wetted            = 653. - 12.*2.4*2 # 2.4 is engine area on one side
    wing.areas.exposed           = 326.5  * Units['meter**2']  
    wing.areas.affected          = .6 * wing.areas.reference
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees
    wing.origin                  = [[14,0,-.8]] # meters
    wing.aerodynamic_center      = [35,0,0] # meters
    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = True
    wing.vortex_lift             = True
    wing.high_mach               = True
    wing.dynamic_pressure_ratio  = 1.0 
    
    # add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'    
    
    wing.aspect_ratio            = 0.74   
    wing.sweeps.quarter_chord    = 60 * Units.deg
    wing.thickness_to_chord      = 0.04
    wing.taper                   = 0.14
    wing.spans.projected         = 6.0   * Units.meter   
    wing.chords.root             = 14.5  * Units.meter
    wing.total_length            = 14.5  * Units.meter
    wing.chords.tip              = 2.7   * Units.meter
    wing.chords.mean_aerodynamic = 8.66  * Units.meter
    wing.areas.reference         = 33.91 * Units['meter**2']  
    wing.areas.wetted            = 76.   * Units['meter**2']  
    wing.areas.exposed           = 38.   * Units['meter**2']  
    wing.areas.affected          = 33.91 * Units['meter**2']  
    wing.twists.root             = 0.0   * Units.degrees
    wing.twists.tip              = 0.0   * Units.degrees  
    wing.origin                  = [[42.,0,1.]] # meters
    wing.aerodynamic_center      = [50,0,0] # meters
    wing.vertical                = True 
    wing.symmetric               = False
    wing.t_tail                  = False
    wing.high_mach               = True     
    wing.dynamic_pressure_ratio  = 1.0
    
    # add to vehicle
    vehicle.append_component(wing)    

    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'
    fuselage.seats_abreast                      = 4
    fuselage.seat_pitch                         = 1     * Units.meter
    fuselage.fineness.nose                      = 4.3   * Units.meter   
    fuselage.fineness.tail                      = 6.4   * Units.meter   
    fuselage.lengths.total                      = 61.66 * Units.meter    
    fuselage.width                              = 2.88  * Units.meter   
    fuselage.heights.maximum                    = 3.32  * Units.meter    
    fuselage.heights.at_quarter_length          = 3.32 * Units.meter   
    fuselage.heights.at_wing_root_quarter_chord = 3.32 * Units.meter   
    fuselage.heights.at_three_quarters_length   = 3.32 * Units.meter   
    fuselage.areas.wetted                       = 447. * Units['meter**2'] 
    fuselage.areas.front_projected              = 11.9 * Units['meter**2'] 
    fuselage.effective_diameter                 = 3.1 * Units.meter    
    fuselage.differential_pressure              = 7.4e4 * Units.pascal    # Maximum differential pressure
    
    # add to vehicle
    vehicle.append_component(fuselage)
    
    # ------------------------------------------------------------------        
    # the nacelle 
    # ------------------------------------------------------------------   
    nacelle                  = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.diameter         = 1.3
    nacelle.tag              = 'nacelle_L1'
    nacelle.origin           = [[36.56, 22, -1.9]] 
    nacelle.length           = 12.0 
    nacelle.inlet_diameter   = 1.1 
    nacelle.areas.wetted     = 30.
    vehicle.append_component(nacelle)       

    nacelle_2               = deepcopy(nacelle)
    nacelle_2.tag           = 'nacelle_2'
    nacelle_2.origin        = [[37.,5.3,-1.3]]
    vehicle.append_component(nacelle_2)     

    nacelle_3               = deepcopy(nacelle)
    nacelle_3.tag           = 'nacelle_3'
    nacelle_3.origin        = [[37.,-5.3,-1.3]]
    vehicle.append_component(nacelle_3)   

    nacelle_4              = deepcopy(nacelle)
    nacelle_4.tag          = 'nacelle_4'
    nacelle_4.origin       = [[37.,-6.,-1.3]]
    vehicle.append_component(nacelle_4)        
         

    # ################################################# Energy Network ##################################################################          
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Turbofjet Network
    #------------------------------------------------------------------------------------------------------------------------------------    
    net                                         = RCAIDE.Energy.Networks.Turbojet_Engine() 

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

    turbojet                                    = RCAIDE.Energy.Converters.Turbojet() 
    turbojet.tag                                = 'turbofan'
    turbojet.origin                             = [[37.,6.,-1.3]] 
    turbojet.engine_length                      = 12.0    
    turbojet.design_altitude                    = 0.0*Units.ft
    turbojet.design_mach_number                 = 0.01
    turbojet.design_thrust                      = 140000. * Units.N  
 
          
    # working fluid
    turbojet.working_fluid = RCAIDE.Attributes.Gases.Air()

    #  Ram 
    ram     = RCAIDE.Energy.Converters.Ram()
    ram.tag = 'ram' 
    turbojet.append(ram)

    # Inlet Nozzle 
    inlet_nozzle = RCAIDE.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag = 'inlet_nozzle' 
    inlet_nozzle.polytropic_efficiency = 0.98
    inlet_nozzle.pressure_ratio        = 1.0 
    turbojet.append(inlet_nozzle)
    
    # Low Pressure Compressor 
    compressor                       = RCAIDE.Energy.Converters.Compressor()    
    compressor.tag                   = 'low_pressure_compressor' 
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 3.1     
    turbojet.append(compressor)
 
    # High Pressure Compressor 
    compressor = RCAIDE.Energy.Converters.Compressor()    
    compressor.tag = 'high_pressure_compressor' 
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 5.0  
    turbojet.append(compressor)
 
    # Low Pressure Turbine 
    turbine = RCAIDE.Energy.Converters.Turbine()   
    turbine.tag='low_pressure_turbine' 
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93      
    turbojet.append(turbine)
     
    # High Pressure Turbine 
    turbine = RCAIDE.Energy.Converters.Turbine()   
    turbine.tag='high_pressure_turbine' 
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93      
    turbojet.append(turbine)
       
    # Combustor 
    combustor = RCAIDE.Energy.Converters.Combustor()   
    combustor.tag = 'combustor' 
    combustor.efficiency                = 0.99   
    combustor.turbine_inlet_temperature = 1450.
    combustor.pressure_ratio            = 1.0
    combustor.fuel_data                 = RCAIDE.Attributes.Propellants.Jet_A()     
    turbojet.append(combustor)
 
    # Core Nozzle 
    nozzle = RCAIDE.Energy.Converters.Supersonic_Nozzle()   
    nozzle.tag = 'core_nozzle' 
    nozzle.polytropic_efficiency = 0.95
    nozzle.pressure_ratio        = 0.99     
    turbojet.append(nozzle) 
 
    # design turbofan
    design_turbojet(turbojet)  
    fuel_line.turbojets.append(turbojet)

    turbojet_2                      = deepcopy(turbojet)
    turbojet_2.tag                  = 'turbojet_2' 
    turbojet_2.origin               = [[37.,5.3,-1.3]] 
    fuel_line.turbojets.append(turbojet_2)      

    turbojet_3                      = deepcopy(turbojet)
    turbojet_3.tag                  = 'turbojet_3' 
    turbojet_3.origin               = [[37.,-5.3,-1.3]] 
    fuel_line.turbojets.append(turbojet_3)    
     

    turbojet_4                      = deepcopy(turbojet)
    turbojet_4.tag                  = 'turbojet_4' 
    turbojet_4.origin               = [[37.,-6.,-1.3]] 
    fuel_line.turbojets.append(turbojet_4)        
 

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
    
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------
    
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'cruise'
    
    configs.append(config)
    
    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------
    
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'takeoff'
    
    config.V2_VS_ratio = 1.21
    config.maximum_lift_coefficient = 2.
    
    configs.append(config)
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'landing'

    config.Vref_VS_ratio = 1.23
    config.maximum_lift_coefficient = 2.
    
    configs.append(config)
    
    return configs

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
    
    # Drag Components
    plot_drag_components(results)
    
    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results)
    
    # Plot Velocities 
    plot_aircraft_velocities(results)      

    return
 

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    
    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
     
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments
    
    # base segment
    base_segment = Segments.Segment()
    
    
    # ------------------------------------------------------------------
    #   First Climb Segment
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1" 
    segment.analyses.extend( analyses.base ) 
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 7. * Units.deg 
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 3.05   * Units.km
    segment.air_speed      = 128.6 * Units['m/s']
    segment.climb_rate     = 20.32 * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------
    #   Second Climb Segment
    # ------------------------------------------------------------------    
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end   = 4.57   * Units.km
    segment.air_speed      = 205.8  * Units['m/s']
    segment.climb_rate     = 10.16  * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------
    #   Third Climb Segment: linear Mach
    # ------------------------------------------------------------------    
    
    segment = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 7.60   * Units.km
    segment.mach_start   = 0.64
    segment.mach_end     = 1.0
    segment.climb_rate   = 5.05  * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: linear Mach
    # ------------------------------------------------------------------    
    
    segment = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_4" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 15.24   * Units.km
    segment.mach_start   = 1.0
    segment.mach_end     = 2.02
    segment.climb_rate   = 5.08  * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    

    # ------------------------------------------------------------------
    #   Fourth Climb Segment
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_5" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 18.288   * Units.km
    segment.mach_number  = 2.02
    segment.climb_rate   = 0.65  * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------    
    #   Cruise Segment
    # ------------------------------------------------------------------    
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.base ) 
    segment.mach       = 2.02
    segment.distance   = 2000.0 * Units.km 
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)        
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------    
    #   First Descent Segment
    # ------------------------------------------------------------------    
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 6.8   * Units.km
    segment.mach_start   = 2.02
    segment.mach_end     = 1.0
    segment.descent_rate = 5.0   * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------    
    #   Second Descent Segment
    # ------------------------------------------------------------------    
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "descent_2" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 3.0   * Units.km
    segment.mach_start   = 1.0
    segment.mach_end     = 0.65
    segment.descent_rate = 5.0   * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------    
    #   Third Descent Segment
    # ------------------------------------------------------------------    

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 130.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']
    segment = analyses.base.energy.networks.turbofan_jet_engine.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------    
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    return mission

def missions_setup(mission):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""

    missions     = RCAIDE.Analyses.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions  

if __name__ == '__main__':  
    main()
    plt.show()