# tut_Concorde.py
 

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE  
from RCAIDE.Framework.Core                                  import Units , Data    
from RCAIDE.Library.Methods.Propulsors.Turbojet_Propulsor   import design_turbojet
from RCAIDE.Library.Methods.Geometry.Planform               import wing_segmented_planform
from RCAIDE.Library.Plots     import *       

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

    # plot vehicle 
    plot_3d_vehicle(vehicle, 
                    min_x_axis_limit            = 0,
                    max_x_axis_limit            = 60,
                    min_y_axis_limit            = -30,
                    max_y_axis_limit            = 30,
                    min_z_axis_limit            = -30,
                    max_z_axis_limit            = 30, 
                    )             
    
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
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)
    
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics                                       = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle                              = vehicle
    aerodynamics.settings.number_of_spanwise_vortices  = 5
    aerodynamics.settings.number_of_chordwise_vortices = 2       
    aerodynamics.settings.model_fuselage               = True
    aerodynamics.settings.drag_coefficient_increment   = 0.0000
    analyses.append(aerodynamics)


    # ------------------------------------------------------------------
    #  Emissions
    emissions = RCAIDE.Framework.Analyses.Emissions.Emission_Index_Correlation_Method()
    emissions.vehicle = vehicle          
    analyses.append(emissions)
  
    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle  = vehicle 
    analyses.append(energy)
    
    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Earth()
    analyses.append(planet)
    
    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   
    
    # done!
    return analyses    

 
def vehicle_setup(): 
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Concorde'
    
    # mass properties
    vehicle.mass_properties.max_takeoff     = 185000.   # kg
    vehicle.mass_properties.operating_empty = 78700.   # kg
    vehicle.mass_properties.takeoff         = 183000.   # kg, adjusted due to significant fuel burn on runway
    vehicle.mass_properties.cargo           = 1000.  * Units.kilogram   
    vehicle.mass_properties.max_zero_fuel   = 92000.
        
    # envelope properties
    vehicle.flight_envelope.ultimate_load   = 3.75
    vehicle.flight_envelope.limit_load      = 2.5
  
    # basic parameters  
    vehicle.reference_area                 = 358.25      
    vehicle.passengers                     = 100
    vehicle.systems.control                = "fully powered" 
    vehicle.systems.accessories            = "sst"
    vehicle.maximum_cross_sectional_area   = 13.9
    vehicle.total_length                   = 61.66
    vehicle.design_mach_number             = 2.02
    vehicle.design_range                   = 4505 * Units.miles
    vehicle.design_cruise_alt              = 60000.0 * Units.ft
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------     
    
    wing = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.aspect_ratio            = 1.83
    wing.sweeps.quarter_chord    = 59.5 * Units.deg
    wing.sweeps.leading_edge     = 66.5 * Units.deg
    wing.thickness_to_chord      = 0.03
    wing.taper                   = 0.
    
    wing.spans.projected           = 25.6    
    
    wing.chords.root               = 33.8
    wing.total_length              = 33.8
    wing.chords.tip                = 1.1
    wing.chords.mean_aerodynamic   = 18.4
    
    wing.areas.reference           = 358.25 
    wing.areas.wetted              = 601.
    wing.areas.exposed             = 326.5
    wing.areas.affected            = .6*wing.areas.reference
    
    wing.twists.root               = 0.0 * Units.degrees
    wing.twists.tip                = 0.0 * Units.degrees
    
    wing.origin                    = [[14,0,-.8]]
    wing.aerodynamic_center        = [35,0,0] 
    
    wing.vertical                  = False
    wing.symmetric                 = True
    wing.high_lift                 = True
    wing.vortex_lift               = True
    wing.high_mach                 = True 
    wing.dynamic_pressure_ratio    = 1.0
     
    wing_airfoil                   = RCAIDE.Library.Components.Airfoils.Airfoil()  
    ospath                         = os.path.abspath(__file__)
    separator                      = os.path.sep
    rel_path                       = os.path.dirname(ospath) + separator   
    wing_airfoil.coordinate_file   = rel_path + 'Airfoils' + separator + 'NACA65_203.txt' 
    wing.append_airfoil(wing_airfoil)  
    
    # set root sweep with inner section
    segment = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'section_1'
    segment.percent_span_location = 0.
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 1
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 67. * Units.deg
    segment.thickness_to_chord    = 0.03
    segment.append_airfoil(wing_airfoil)
    wing.Segments.append(segment)
    
    # set section 2 start point
    segment = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'section_2'
    segment.percent_span_location = (6.15 * 2) /wing.spans.projected
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 13.8/wing.chords.root
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 48. * Units.deg
    segment.thickness_to_chord    = 0.03
    segment.append_airfoil(wing_airfoil)
    wing.Segments.append(segment)
    
    
    # set section 3 start point
    segment = RCAIDE.Library.Components.Wings.Segment() 
    segment.tag                   = 'section_3'
    segment.percent_span_location = (12.1 *2) /wing.spans.projected
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 4.4/wing.chords.root
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 71. * Units.deg 
    segment.thickness_to_chord    = 0.03
    segment.append_airfoil(wing_airfoil)
    wing.Segments.append(segment)  
    
    # set tip
    segment = RCAIDE.Library.Components.Wings.Segment() 
    segment.tag                   = 'tip'
    segment.percent_span_location = 1.
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 1.1/wing.chords.root
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 0.
    segment.thickness_to_chord    = 0.03
    segment.append_airfoil(wing_airfoil)
    wing.Segments.append(segment)      
    
    # Fill out more segment properties automatically
    wing = wing_segmented_planform(wing)        
    
    
    # add to vehicle
    vehicle.append_component(wing)
    
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'    
    
    wing.aspect_ratio            = 0.74     
    wing.sweeps.quarter_chord    = 60 * Units.deg
    wing.thickness_to_chord      = 0.04
    wing.taper                   = 0.14 
    wing.spans.projected         = 6.0     
    wing.chords.root             = 14.5
    wing.total_length            = 14.5
    wing.chords.tip              = 2.7
    wing.chords.mean_aerodynamic = 8.66 
    wing.areas.reference         = 33.91     
    wing.areas.wetted            = 76. 
    wing.areas.exposed           = 38.
    wing.areas.affected          = 33.91 
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees   
    wing.origin                  = [[42.,0,1.]]
    wing.aerodynamic_center      = [50,0,0]     
    wing.vertical                = True 
    wing.symmetric               = False
    wing.t_tail                  = False
    wing.high_mach               = True     
    
    wing.dynamic_pressure_ratio  = 1.0
    
    tail_airfoil = RCAIDE.Library.Components.Airfoils.Airfoil() 
    tail_airfoil.coordinate_file = rel_path + 'Airfoils' + separator + 'supersonic_tail.txt' 
    
    wing.append_airfoil(tail_airfoil)  

    # set root sweep with inner section
    segment = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'section_1'
    segment.percent_span_location = 0.0
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 14.5/14.5
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 63. * Units.deg
    segment.thickness_to_chord    = 0.04
    segment.append_airfoil(tail_airfoil)
    wing.Segments.append(segment)
    
    # set mid section start point
    segment = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'section_2'
    segment.percent_span_location = 2.4/(6.0) + wing.Segments['section_1'].percent_span_location
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 7.5/14.5
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 40. * Units.deg
    segment.thickness_to_chord    = 0.04
    segment.append_airfoil(tail_airfoil)
    wing.Segments.append(segment)
    
    # set tip
    segment = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'tip'
    segment.percent_span_location = 1.
    segment.twist                 = 0. * Units.deg
    segment.root_chord_percent    = 2.7/14.5
    segment.dihedral_outboard     = 0.
    segment.sweeps.quarter_chord  = 0.
    segment.thickness_to_chord    = 0.04
    segment.append_airfoil(tail_airfoil)
    wing.Segments.append(segment)    
    
    # Fill out more segment properties automatically
    wing = wing_segmented_planform(wing)        
    
    # add to vehicle
    vehicle.append_component(wing)    


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage                                        = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.seats_abreast                          = 4
    fuselage.seat_pitch                             = 38. * Units.inches 
    fuselage.fineness.nose                          = 4.3
    fuselage.fineness.tail                          = 6.4 
    fuselage.lengths.total                          = 61.66   
    fuselage.width                                  = 2.88 
    fuselage.heights.maximum                        = 3.32    
    fuselage.heights.maximum                        = 3.32    
    fuselage.heights.at_quarter_length              = 3.32    
    fuselage.heights.at_wing_root_quarter_chord     = 3.32    
    fuselage.heights.at_three_quarters_length       = 3.32    
    fuselage.areas.wetted                           = 442.
    fuselage.areas.front_projected                  = 11.9 
    fuselage.effective_diameter                     = 3.1 
    fuselage.differential_pressure                  = 7.4e4 * Units.pascal    # Maximum differential pressure 
    
    fuselage.OpenVSP_values = Data() # VSP uses degrees directly
    
    fuselage.OpenVSP_values.nose = Data()
    fuselage.OpenVSP_values.nose.top = Data()
    fuselage.OpenVSP_values.nose.side = Data()
    fuselage.OpenVSP_values.nose.top.angle = 20.0
    fuselage.OpenVSP_values.nose.top.strength = 0.75
    fuselage.OpenVSP_values.nose.side.angle = 20.0
    fuselage.OpenVSP_values.nose.side.strength = 0.75  
    fuselage.OpenVSP_values.nose.TB_Sym = True
    fuselage.OpenVSP_values.nose.z_pos = -.01
    
    fuselage.OpenVSP_values.tail = Data()
    fuselage.OpenVSP_values.tail.top = Data()
    fuselage.OpenVSP_values.tail.side = Data()    
    fuselage.OpenVSP_values.tail.bottom = Data()
    fuselage.OpenVSP_values.tail.top.angle = 0.0
    fuselage.OpenVSP_values.tail.top.strength = 0.0 
    
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_0'    
    segment.percent_x_location                  = 0.0000
    segment.percent_z_location                  =  -0.61 /fuselage.lengths.total 
    segment.height                              = 0.00100 
    segment.width                               = 0.00100  
    fuselage.Segments.append(segment)   
     
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'   
    segment.percent_x_location                  = 3.02870/fuselage.lengths.total   
    segment.percent_z_location                  = -0.3583/fuselage.lengths.total     
    segment.height                              = 1.4502  
    segment.width                               = 1.567  
    fuselage.Segments.append(segment)
    

    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'   
    segment.percent_x_location                  =   5.7742/fuselage.lengths.total   
    segment.percent_z_location                  =  -0.1500/fuselage.lengths.total    
    segment.height                              = 2.356  
    segment.width                               = 2.3429  
    fuselage.Segments.append(segment)
    

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  =  9.0791/fuselage.lengths.total    
    segment.percent_z_location                  = 0  
    segment.height                              = 3.0581  
    segment.width                               = 2.741  
    fuselage.Segments.append(segment)
     
    
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_4'    
    segment.percent_x_location                  = 12.384/fuselage.lengths.total  
    segment.percent_z_location                  = 0  
    segment.height                              = 3.3200 
    segment.width                               = 2.880  
    fuselage.Segments.append(segment)
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 43.228 /fuselage.lengths.total    
    segment.percent_z_location                  = 0 
    segment.height                              = 3.3200   
    segment.width                               = 2.8800  
    fuselage.Segments.append(segment)

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  =  47.5354/fuselage.lengths.total     
    segment.percent_z_location                  =  0.100/fuselage.lengths.total     
    segment.height                              = 2.952  
    segment.width                               = 2.8800  
    fuselage.Segments.append(segment)   

                 
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 1  
    segment.percent_z_location                  = 1.2332/fuselage.lengths.total   
    segment.height                              = 0.00100  
    segment.width                               = 0.00100  
    fuselage.Segments.append(segment)
    
    
    vehicle.append_component(fuselage)

    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################## Energy Network ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    #initialize the fuel network
    net                                            = RCAIDE.Framework.Networks.Fuel() 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                     = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line() 
    fuel_line.identical_propulsors                = False # for regression 
    

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Right Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    outer_right_turbojet                          = RCAIDE.Library.Components.Propulsors.Turbojet()  
    outer_right_turbojet.tag                      = 'outer_right_turbojet'   
    outer_right_turbojet.active_fuel_tanks        = ['tank_6_and_7','tank_5A_and_7A','tank_2_and_3','tank_11']    
    outer_right_turbojet.engine_length            = 4.039
    outer_right_turbojet.nacelle_diameter         = 1.3
    outer_right_turbojet.inlet_diameter           = 1.212 
    outer_right_turbojet.areas.wetted             = 30
    outer_right_turbojet.design_altitude          = 60000.0*Units.ft
    outer_right_turbojet.design_mach_number       = 2.02
    outer_right_turbojet.design_thrust            = 10000. * Units.lbf  
    outer_right_turbojet.origin                   = [[37.,5.5,-1.6]] 
    outer_right_turbojet.working_fluid            = RCAIDE.Library.Attributes.Gases.Air()
    
    # Ram  
    ram                                           = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                       = 'ram' 
    outer_right_turbojet.ram                      = ram 
         
    # Inlet Nozzle         
    inlet_nozzle                                  = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                              = 'inlet_nozzle' 
    inlet_nozzle.polytropic_efficiency            = 1.0
    inlet_nozzle.pressure_ratio                   = 1.0
    inlet_nozzle.pressure_recovery                = 0.94 
    outer_right_turbojet.inlet_nozzle             = inlet_nozzle    
          
    #  Low Pressure Compressor      
    lp_compressor                                 = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    lp_compressor.tag                             = 'low_pressure_compressor' 
    lp_compressor.polytropic_efficiency           = 0.88
    lp_compressor.pressure_ratio                  = 3.1     
    outer_right_turbojet.low_pressure_compressor  = lp_compressor         
        
    # High Pressure Compressor        
    hp_compressor                                 = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    hp_compressor.tag                             = 'high_pressure_compressor' 
    hp_compressor.polytropic_efficiency           = 0.88
    hp_compressor.pressure_ratio                  = 5.0  
    outer_right_turbojet.high_pressure_compressor = hp_compressor
 
    # Low Pressure Turbine 
    lp_turbine                                    = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    lp_turbine.tag                                ='low_pressure_turbine' 
    lp_turbine.mechanical_efficiency              = 0.99
    lp_turbine.polytropic_efficiency              = 0.89 
    outer_right_turbojet.low_pressure_turbine     = lp_turbine      
             
    # High Pressure Turbine         
    hp_turbine                                    = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    hp_turbine.tag                                ='high_pressure_turbine' 
    hp_turbine.mechanical_efficiency              = 0.99
    hp_turbine.polytropic_efficiency              = 0.87 
    outer_right_turbojet.high_pressure_turbine    = hp_turbine   
          
    # Combustor   
    combustor                                     = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                 = 'combustor' 
    combustor.efficiency                          = 0.94
    combustor.alphac                              = 1.0     
    combustor.turbine_inlet_temperature           = 1440.
    combustor.pressure_ratio                      = 0.92
    combustor.fuel_data                           = RCAIDE.Library.Attributes.Propellants.Jet_A()     
    outer_right_turbojet.combustor                = combustor
     
    #  Afterburner  
    afterburner                                   = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    afterburner.tag                               = 'afterburner' 
    afterburner.efficiency                        = 0.9
    afterburner.alphac                            = 1.0     
    afterburner.turbine_inlet_temperature         = 1500
    afterburner.pressure_ratio                    = 1.0
    afterburner.fuel_data                         = RCAIDE.Library.Attributes.Propellants.Jet_A()     
    outer_right_turbojet.afterburner              = afterburner   
 
    # Core Nozzle 
    nozzle                                        = RCAIDE.Library.Components.Propulsors.Converters.Supersonic_Nozzle()   
    nozzle.tag                                    = 'core_nozzle' 
    nozzle.pressure_recovery                      = 0.95
    nozzle.pressure_ratio                         = 1.    
    outer_right_turbojet.core_nozzle              = nozzle
    
    # design turbojet 
    design_turbojet(outer_right_turbojet) 

    nacelle                                     = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    nacelle.diameter                            = 1.3
    nacelle.tag                                 = 'nacelle_1'
    nacelle.origin                              = [[37.,5.5,-1.6]] 
    nacelle.length                              = 10
    nacelle.inlet_diameter                      = 1.1 
    nacelle.areas.wetted                        = 30.
    
    nac_segment                                 = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                             = 'segment_1' 
    nac_segment.orientation_euler_angles        = [0., -45*Units.degrees,0.]     
    nac_segment.percent_x_location              = 0.0  
    nac_segment.height                          = 2.12
    nac_segment.width                           = 1.5
    nac_segment.curvature                       = 5
    nacelle.append_segment(nac_segment)         

    nac_segment                                 = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                             = 'segment_2'
    nac_segment.percent_x_location              = 1.0
    nac_segment.height                          = 1.5
    nac_segment.width                           = 1.5
    nac_segment.curvature                       = 5
    nacelle.append_segment(nac_segment)      
    outer_right_turbojet.nacelle = nacelle  
    net.propulsors.append(outer_right_turbojet) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Right Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    inner_right_turbojet                     = deepcopy(outer_right_turbojet) 
    inner_right_turbojet.tag                 = 'inner_right_turbojet'      
    inner_right_turbojet.origin              = [[37.,4,-1.6]]   
    inner_right_turbojet.active_fuel_tanks   = ['tank_6_and_7','tank_5A_and_7A','tank_2_and_3','tank_11']  
    nacelle_2                                = deepcopy(nacelle)
    nacelle_2.tag                            = 'nacelle_2'
    nacelle_2.origin                         = [[37.,4,-1.6]]
    inner_right_turbojet.nacelle = nacelle_2 
    net.propulsors.append(inner_right_turbojet) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Right Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------    
    inner_left_turbojet                     = deepcopy(outer_right_turbojet)    
    inner_left_turbojet.tag                 = 'inner_left_turbojet'  
    inner_left_turbojet.origin              = [[37.,-4,-1.6]]  
    inner_left_turbojet.active_fuel_tanks   = ['tank_9','tank_10','tank_1_and_4','tank_5_and_8'] 
    nacelle_3                               = deepcopy(nacelle)
    nacelle_3.tag                           = 'nacelle_3'
    nacelle_3.origin                        = [[37.,-4,-1.6]]
    inner_left_turbojet.nacelle = nacelle_3 
    net.propulsors.append(inner_left_turbojet) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Left Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------    
    outer_left_turbojet                     = deepcopy(outer_right_turbojet)
    outer_left_turbojet.tag                 = 'outer_left_turbojet'     
    outer_left_turbojet.active_fuel_tanks   = ['tank_9','tank_10','tank_1_and_4','tank_5_and_8']   
    outer_left_turbojet.origin              = [[37.,-5.5,-1.6]]   
    nacelle_4                               = deepcopy(nacelle)
    nacelle_4.tag                           = 'nacelle_4'
    nacelle_4.origin                        = [[37.,-5.5,-1.6]]
    outer_left_turbojet.nacelle = nacelle_4
    net.propulsors.append(outer_left_turbojet) 
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------   
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_9'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[26.5,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 11096
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                            = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_10'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[28.7,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 11943
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                            = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_1_and_4'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[31.0,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 4198+4198
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                            = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_5_and_8'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[32.9,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 7200+12838
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                              = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_6_and_7'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[37.4,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 11587+7405
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                                 = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_5A_and_7A'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[40.2,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 2225+2225
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                                 = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_2_and_3'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[40.2,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 4570+4570
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                                 = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank)  
 
    fuel_tank = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_11'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[49.8,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 10415
    fuel_tank.fuel_selector_ratio                  = 1/8
    fuel_tank.fuel                            = RCAIDE.Library.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank)
     
    #------------------------------------------------------------------------------------------------------------------------------------   
    # Assign propulsors to fuel line to network      
    fuel_line.assigned_propulsors =  [[outer_left_turbojet.tag,inner_left_turbojet.tag, outer_right_turbojet.tag, inner_right_turbojet.tag]]    
    
     # Append fuel line to network      
    net.fuel_lines.append(fuel_line)    
  
    #------------------------------------------------------------------------------------------------------------------------------------          
    # Append energy network to aircraft 
    vehicle.append_energy_network(net)     

    return vehicle


# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs         = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config     = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------ 
    config     = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cruise' 
    configs.append(config)
    
    # ------------------------------------------------------------------
    #   Afterburner Climb Configuration
    # ------------------------------------------------------------------ 
    config                                      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                                  = 'climb' 
    for propulsor in config.networks.fuel.propulsors:
        propulsor.afterburner_active = True 
    configs.append(config)    
    
    
    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------  
    config                                      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                                  = 'takeoff'  
    config.maximum_lift_coefficient             = 2.  
    for propulsor in config.networks.fuel.propulsors:
        propulsor.afterburner_active = True 
    configs.append(config) 
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config                                = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                            = 'landing' 
    config.wings['main_wing'].flaps_angle = 0. * Units.deg
    config.wings['main_wing'].slats_angle = 0. * Units.deg  
    config.maximum_lift_coefficient       = 2. 
    
    configs.append(config)
    
    # done!
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
    
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
     
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments 
    base_segment = Segments.Segment()
    
    # ------------------------------------------------------------------
    #   First Climb Segment: constant Mach, constant segment angle 
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1" 
    segment.analyses.extend( analyses.climb ) 
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 4000. * Units.ft
    segment.air_speed      = 250.  * Units.kts
    segment.climb_rate     = 3000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------
    #   Second Climb Segment: constant Speed, constant segment angle 
    # ------------------------------------------------------------------    
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end = 8000. * Units.ft
    segment.air_speed    = 250.  * Units.kts
    segment.climb_rate   = 2000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Second Climb Segment: constant Speed, constant segment angle 
    # ------------------------------------------------------------------    
    
    segment = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end        = 33000. * Units.ft
    segment.mach_number_start   = .45
    segment.mach_number_end     = 0.95
    segment.climb_rate          = 3000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)    

    # ------------------------------------------------------------------
    #   Third Climb Segment: linear Mach, constant segment angle 
    # ------------------------------------------------------------------    
      
    segment = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.climb ) 
    segment.altitude_end        = 34000. * Units.ft
    segment.mach_number_start   = 0.95
    segment.mach_number_end     = 1.1
    segment.climb_rate          = 2000.  * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment) 

    # ------------------------------------------------------------------
    #   Third Climb Segment: linear Mach, constant segment angle 
    # ------------------------------------------------------------------   
    segment     = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_4" 
    segment.analyses.extend( analyses.climb ) 
    segment.altitude_end        = 40000. * Units.ft
    segment.mach_number_start   = 1.1
    segment.mach_number_end     = 1.7
    segment.climb_rate          = 1750.  * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: linear Mach, constant segment angle 
    # ------------------------------------------------------------------  
    segment = Segments.Climb.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "climb_5" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end        = 50000. * Units.ft
    segment.mach_number_start   = 1.7
    segment.mach_number_end     = 2.02
    segment.climb_rate          = 750.  * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)     
    

    # ------------------------------------------------------------------
    #   Fourth Climb Segment: linear Mach, constant segment angle 
    # ------------------------------------------------------------------  
    segment = Segments.Climb.Constant_Mach_Constant_Rate(base_segment)
    segment.tag = "climbing_cruise" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude_end = 56500. * Units.ft
    segment.mach_number  = 2.02
    segment.climb_rate   = 50.  * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------    
    #   Cruise Segment: constant speed 
    # ------------------------------------------------------------------    
    segment     = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "level_cruise" 
    segment.analyses.extend( analyses.cruise ) 
    segment.mach_number                          = 2.02
    segment.distance                             = 1. * Units.nmi
    segment.state.numerics.number_control_points = 4  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   First Descent Segment: decceleration
    # ------------------------------------------------------------------    
    segment     = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag = "decel_1" 
    segment.analyses.extend( analyses.cruise )
    segment.acceleration      = -.5  * Units['m/s/s']
    segment.air_speed_end     = 1.5*573.  * Units.kts  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   First Descent Segment
    # ------------------------------------------------------------------  
    segment     = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise )
    segment.altitude_end      = 41000. * Units.ft
    segment.mach_number_end   = 1.3
    segment.descent_rate      = 2000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #   First Descent Segment: decceleration
    # ------------------------------------------------------------------   
    segment     = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag = "decel_2" 
    segment.analyses.extend( analyses.cruise )
    segment.acceleration      = -.5  * Units['m/s/s']
    segment.air_speed_end     = 0.95*573.  * Units.kts  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #   First Descent Segment
    # ------------------------------------------------------------------    
      
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "descent_2" 
    segment.analyses.extend( analyses.cruise )
    segment.altitude_end      = 10000. * Units.ft
    segment.mach_number_end   = 250./638. 
    segment.descent_rate      = 2000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #   First Descent Segment
    # ------------------------------------------------------------------     
    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3" 
    segment.analyses.extend( analyses.cruise )
    segment.altitude_end = 0. * Units.ft
    segment.air_speed    = 250. * Units.kts
    segment.descent_rate = 1000. * Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['inner_right_turbojet','outer_right_turbojet','outer_left_turbojet','inner_left_turbojet']]
    segment.assigned_control_variables.throttle.initial_values       = [[0.5]] 
    segment.assigned_control_variables.body_angle.active             = True               
    segment.assigned_control_variables.body_angle.initial_values     = [[3*Units.degrees]]    
    
    mission.append_segment(segment)      
    
    # ------------------------------------------------------------------    
    #   Mission definition complete    
    # ------------------------------------------------------------------ 
    return mission

def missions_setup(mission):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""

    missions     = RCAIDE.Framework.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions  

if __name__ == '__main__':  
    main()
    plt.show()