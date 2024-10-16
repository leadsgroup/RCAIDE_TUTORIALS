# tut_solar_UAV.py
# 
# Created: Dec 2021, M.Clarke

#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
import RCAIDE
from RCAIDE.Library.Plots                                    import *    
from RCAIDE.Framework.Core                                   import Units 
from RCAIDE.Framework.Networks.Solar                         import Solar
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor      import design_propeller
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common  import initialize_from_mass

import matplotlib.pyplot as plt
import numpy as np 
import time

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
    mission  = mission_setup(analyses) 

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
    results = missions.base.evaluate()  
       
    # plot results    
    plot_mission(results)
    
    return 

# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Solar'
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    
    # mass properties
    vehicle.mass_properties.takeoff         = 250. * Units.kg
    vehicle.mass_properties.operating_empty = 250. * Units.kg
    vehicle.mass_properties.max_takeoff     = 250. * Units.kg 
    
    # basic parameters
    vehicle.reference_area                    = 80.       
    vehicle.envelope.ultimate_load            = 2.0
    vehicle.envelope.limit_load               = 1.5
    vehicle.envelope.maximum_dynamic_pressure = 0.5*1.225*(40.**2.) #Max q

    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------   

    wing = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.areas.reference         = vehicle.reference_area
    wing.spans.projected         = 40.0 * Units.meter
    wing.aspect_ratio            = (wing.spans.projected**2)/wing.areas.reference 
    wing.sweeps.quarter_chord    = 0.0 * Units.deg
    wing.symmetric               = True
    wing.thickness_to_chord      = 0.12
    wing.taper                   = 1.0
    wing.vertical                = False
    wing.high_lift               = True 
    wing.dynamic_pressure_ratio  = 1.0
    wing.chords.mean_aerodynamic = wing.areas.reference/wing.spans.projected
    wing.chords.root             = wing.areas.reference/wing.spans.projected
    wing.chords.tip              = wing.areas.reference/wing.spans.projected
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees
    wing.highlift                = False  
    wing.vertical                = False 
    wing.number_ribs             = 26.
    wing.number_end_ribs         = 2.
    wing.transition_x_upper      = 0.6
    wing.transition_x_lower      = 1.0
    wing.origin                  = [[3.0,0.0,0.0]] # meters
    wing.aerodynamic_center      = [1.0,0.0,0.0] # meters
    
    # add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------        
    
    wing = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'
    
    wing.aspect_ratio         = 20. 
    wing.sweeps.quarter_chord = 0 * Units.deg
    wing.thickness_to_chord   = 0.12
    wing.taper                = 1.0
    wing.areas.reference      = vehicle.reference_area * .15
    wing.areas.wetted         = 2.0 * wing.areas.reference
    wing.areas.exposed        = 0.8 * wing.areas.wetted
    wing.areas.affected       = 0.6 * wing.areas.wetted       
    wing.spans.projected      = np.sqrt(wing.aspect_ratio*wing.areas.reference)
    wing.twists.root          = 0.0 * Units.degrees
    wing.twists.tip           = 0.0 * Units.degrees      
    
    wing.vertical                = False 
    wing.symmetric               = True
    wing.dynamic_pressure_ratio  = 0.9      
    wing.number_ribs             = 5.0
    wing.chords.root             = wing.areas.reference/wing.spans.projected
    wing.chords.tip              = wing.areas.reference/wing.spans.projected
    wing.chords.mean_aerodynamic = wing.areas.reference/wing.spans.projected  
    wing.origin                  = [[10.,0.0,0.0]] # meters
    wing.aerodynamic_center      = [0.5,0.0,0.0] # meters
  
    # add to vehicle
    vehicle.append_component(wing)    
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'    
    
    
    wing.aspect_ratio         = 20.       
    wing.sweeps.quarter_chord = 0 * Units.deg
    wing.thickness_to_chord   = 0.12
    wing.taper                = 1.0
    wing.areas.reference      = vehicle.reference_area * 0.1
    wing.spans.projected      = np.sqrt(wing.aspect_ratio*wing.areas.reference)

    wing.chords.root             = wing.areas.reference/wing.spans.projected
    wing.chords.tip              = wing.areas.reference/wing.spans.projected
    wing.chords.mean_aerodynamic = wing.areas.reference/wing.spans.projected 
    wing.areas.wetted            = 2.0 * wing.areas.reference
    wing.areas.exposed           = 0.8 * wing.areas.wetted
    wing.areas.affected          = 0.6 * wing.areas.wetted    
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [[10.,0.0,0.0]] # meters
    wing.aerodynamic_center      = [0.5,0.0,0.0] # meters
    wing.symmetric               = True          
    wing.vertical                = True 
    wing.t_tail                  = False
    wing.dynamic_pressure_ratio  = 1.0
    wing.number_ribs             = 5.
  
    # add to vehicle
    vehicle.append_component(wing)  
    
    
    # ------------------------------------------------------------------
    #   Nacelle  
    # ------------------------------------------------------------------
    nacelle              = RCAIDE.Library.Components.Nacelles.Nacelle()
    nacelle.diameter     = 0.2 * Units.meters
    nacelle.length       = 0.01 * Units.meters
    nacelle.tag          = 'nacelle' 
    nacelle.areas.wetted =  nacelle.length *(2*np.pi*nacelle.diameter/2.)
    vehicle.append_component(nacelle) 
        
    
    #------------------------------------------------------------------
    # Propulsor
    #------------------------------------------------------------------
    
    # build network
    net = Solar()  

    #------------------------------------------------------------------------------------------------------------------------------------      
    # Sun 
    #------------------------------------------------------------------------------------------------------------------------------------  
    sun                                            = RCAIDE.Energy.Processes.Solar_Radiation()
    net.solar_flux                                 = sun 
    net.maximum_power_point_tracking_efficency     = 0.95 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                        = RCAIDE.Library.Components.Energy.Distribution.Electrical_Bus()
    bus.fixed_voltage          = True 
    bus.voltage                = 120.0

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Solar Panel 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    panel                      = RCAIDE.Library.Components.Propulsors.Converters.Solar_Panel()
    panel.area                 = vehicle.reference_area * 0.9
    panel.efficiency           = 0.25
    panel.mass_properties.mass = panel.area*(0.60 * Units.kg)
    net.solar_panel            = panel


    #------------------------------------------------------------------------------------------------------------------------------------  
    # Electronic Speed Controller    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    esc            = RCAIDE.Library.Components.Propulsors.Modulators.Electronic_Speed_Controller()
    esc.efficiency = 0.95 
    bus.electronic_speed_controllers.append(esc)  


    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat = RCAIDE.Library.Components.Energy.Batteries.Lithium_Ion_Generic() 
    bat.mass_properties.mass = 95.0 * Units.kg
    bat.specific_energy      = 800. * Units.Wh/Units.kg
    bat.pack.nominal_voltage = 130.0
    initialize_from_mass(bat) 
    bus.batteries.append(bat)                                

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propeller    
    #------------------------------------------------------------------------------------------------------------------------------------    
    propeller                                   = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
    propeller.number_of_blades                  = 2.0
    propeller.tip_radius                        = 4.25 * Units.meters
    propeller.hub_radius                        = 0.05 * Units.meters
    propeller.cruise.design_freestream_velocity = 40.0 * Units['m/s']# freestream
    propeller.cruise.design_angular_velocity    = 150. * Units['rpm']
    propeller.cruise.design_Cl                  = 0.7
    propeller.cruise.design_altitude            = 15.0 * Units.km
    propeller.cruise.design_power               = None
    propeller.cruise.design_thrust              = 120.
    propeller                                   = design_propeller(propeller) 
    bus.rotors.append(propeller)   

    #------------------------------------------------------------------------------------------------------------------------------------           
    # Motor 
    #------------------------------------------------------------------------------------------------------------------------------------   
    motor                      = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor() 
    motor.resistance           = 0.006
    motor.no_load_current      = 2.5  * Units.ampere
    motor.speed_constant       = 30. * Units['rpm'] # RPM/volt converted to (rad/s)/volt    
    motor.rotor_radius         = propeller.tip_radius
    motor.rotor_CP             = propeller.cruise.design_power_coefficient
    motor.gear_ratio           = 12. # Gear ratio
    motor.gearbox_efficiency   = .98 # Gear box efficiency
    motor.expected_current     = 60. # Expected current
    motor.mass_properties.mass = 2.0  * Units.kg 
    bus.motors.append(motor)

    # append bus   
    net.busses.append(bus)

    #------------------------------------------------------------------------------------------------------------------------------------           
    # Payload 
    #------------------------------------------------------------------------------------------------------------------------------------
    payload                      = RCAIDE.Library.Components.Payloads.Payload()
    payload.power_draw           = 50. * Units.watts 
    payload.mass_properties.mass = 5.0 * Units.kg
    net.payload                  = payload

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------  
    avionics            = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw = 50. * Units.watts
    net.avionics        = avionics     

    # add the solar network to the vehicle
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
    base_config.networks.solar.busses.bus.active_propulsor_groups = ['propulsor']     
    base_config.tag = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------
    
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cruise' 
    configs.append(config)
    
    return configs


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
    weights = RCAIDE.Framework.Analyses.Weights.Weights_UAV()
    weights.settings.empty = RCAIDE.Library.Methods.Weights.Correlation_Buildups.UAV.empty
    weights.vehicle = vehicle
    analyses.append(weights)
    
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
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
    mission.tag = 'The Test Mission'

    mission.atmosphere  = RCAIDE.Library.Attributes.Atmospheres.Earth.US_Standard_1976()
    mission.planet      = RCAIDE.Library.Attributes.Planets.Earth()
    
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments
    
    # base segment
    base_segment = Segments.Segment()    
    
    # ------------------------------------------------------------------    
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------    
    
    segment = RCAIDE.Framework.Mission.Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "cruise1"
    
    # connect vehicle configuration
    segment.analyses.extend(analyses.cruise)
    
    # segment attributes     
    segment.start_time                      = time.strptime("Tue, Jun 21 11:30:00  2022", "%a, %b %d %H:%M:%S %Y",)
    segment.altitude                        = 15.0  * Units.km 
    segment.mach_number                     = 0.12
    segment.distance                        = 3050.0 * Units.km
    segment.initial_battery_state_of_charge = 0.3  
    segment.latitude                        = 37.4300   # this defaults to degrees (do not use Units.degrees)
    segment.longitude                       = -122.1700 # this defaults to degrees      
    mission.append_segment(segment)    

    # ------------------------------------------------------------------    
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    return mission

def missions_setup(base_mission):

    # the mission container
    missions = RCAIDE.Framework.Mission.Missions()
    
    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------
    
    missions.base = base_mission
    
    # done!
    return missions  

# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------
def plot_mission(results):

    # Plot Flight Conditions 
    plot_flight_conditions(results) 
    
    # Plot Solar Conditions 
    plot_solar_network_conditions(results)
    
    plot_aerodynamic_forces(results)

    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
    
    # Drag Components
    plot_drag_components(results)    

    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results)

    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results) 
    
    return 

if __name__ == '__main__':
    main()
    plt.show()