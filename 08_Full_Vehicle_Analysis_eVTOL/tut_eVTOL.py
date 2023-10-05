''' 
# Vehicle.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Core import Units, Data   
from RCAIDE.Energy.Networks.All_Electric                                      import All_Electric 
from RCAIDE.Methods.Performance.estimate_cruise_drag                          import estimate_cruise_drag
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from RCAIDE.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration 
from RCAIDE.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from RCAIDE.Methods.Propulsion                                                import size_optimal_motor
from RCAIDE.Methods.Propulsion                                                import design_propeller ,design_lift_rotor 
from RCAIDE.Methods.Weights.Buildups.eVTOL                                    import compute_weight
from RCAIDE.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                         import wing_segmented_planform 
from RCAIDE.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
from RCAIDE.Visualization                 import *       
 
import os
import numpy as np 
from copy import deepcopy  
import matplotlib.pyplot as plt  

# ----------------------------------------------------------------------------------------------------------------------
#  REGRESSION
# ----------------------------------------------------------------------------------------------------------------------  
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
    plot_results(results) 
           
        
    return  


# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup() :
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------  
     
         
    vehicle               = RCAIDE.Vehicle()
    vehicle.tag           = 'Lift_Cruise'
    vehicle.configuration = 'eVTOL'
     
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties #####################################################  
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # mass properties 
    vehicle.mass_properties.max_takeoff       = 2700 
    vehicle.mass_properties.takeoff           = vehicle.mass_properties.max_takeoff
    vehicle.mass_properties.operating_empty   = vehicle.mass_properties.max_takeoff
    vehicle.envelope.ultimate_load            = 5.7   
    vehicle.envelope.limit_load               = 3.  
    vehicle.passengers                        = 6 
        
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Wings ################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------ 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Main Wing
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                          = RCAIDE.Components.Wings.Main_Wing()
    wing.tag                      = 'main_wing'  
    wing.aspect_ratio             = 8.95198  # will  be overwritten
    wing.sweeps.quarter_chord     = 0.0  
    wing.thickness_to_chord       = 0.14 
    wing.taper                    = 0.292
    wing.spans.projected          = 11.82855
    wing.chords.root              = 1.75
    wing.total_length             = 1.75
    wing.chords.tip               = 1.0
    wing.chords.mean_aerodynamic  = 0.8
    wing.dihedral                 = 0.0  
    wing.areas.reference          = 15.629
    wing.twists.root              = 4. * Units.degrees
    wing.twists.tip               = 0. 
    wing.origin                   = [[1.5, 0., 0.991]]
    wing.aerodynamic_center       = [ 1.567, 0., 0.991]    
    wing.winglet_fraction         = 0.0  
    wing.symmetric                = True
    wing.vertical                 = False
    
    ospath                        = os.path.abspath(__file__)
    separator                     = os.path.sep
    rel_path                      = ospath.split( 'tut_eVTOL.py')[0]  +  '..' + separator   
    airfoil                       = RCAIDE.Components.Airfoils.Airfoil()
    airfoil.coordinate_file       = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'
    
    # Segment                                  
    segment                       = RCAIDE.Components.Wings.Segment()
    segment.tag                   = 'Section_1'   
    segment.percent_span_location = 0.0
    segment.twist                 = 4. * Units.degrees 
    segment.root_chord_percent    = 1. 
    segment.dihedral_outboard     = 8 * Units.degrees
    segment.sweeps.quarter_chord  = 0.9  * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)               
    
    # Segment                                   
    segment                       = RCAIDE.Components.Wings.Segment()
    segment.tag                   = 'Section_2'    
    segment.percent_span_location = 3.5/wing.spans.projected
    segment.twist                 = 3. * Units.degrees 
    segment.root_chord_percent    = 1.4000/1.7500
    segment.dihedral_outboard     = 0.0 * Units.degrees
    segment.sweeps.quarter_chord  = 1.27273 * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)               
     
    # Segment                                  
    segment                       = RCAIDE.Components.Wings.Segment()
    segment.tag                   = 'Section_3'   
    segment.percent_span_location = 11.3/wing.spans.projected 
    segment.twist                 = 2.0 * Units.degrees 
    segment.root_chord_percent    = 1.000/1.7500
    segment.dihedral_outboard     = 35.000* Units.degrees 
    segment.sweeps.quarter_chord  = 45.000* Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)     
    
    # Segment                                  
    segment                       = RCAIDE.Components.Wings.Segment()
    segment.tag                   = 'Section_4'   
    segment.percent_span_location = 11.6/wing.spans.projected 
    segment.twist                 = 0.0 * Units.degrees 
    segment.root_chord_percent    = 0.9/1.7500
    segment.dihedral_outboard     = 60. * Units.degrees 
    segment.sweeps.quarter_chord  = 70.0 * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)  
    
    # Segment                                  
    segment                       = RCAIDE.Components.Wings.Segment()
    segment.tag                   = 'Section_5'   
    segment.percent_span_location = 1.0
    segment.twist                 = 0.0 * Units.degrees 
    segment.root_chord_percent    = 0.35/1.7500
    segment.dihedral_outboard     = 0  * Units.degrees 
    segment.sweeps.quarter_chord  = 0  * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)                 
    
    
    # compute reference properties 
    wing_segmented_planform(wing, overwrite_reference = True ) 
    wing = segment_properties(wing)
    vehicle.reference_area        = wing.areas.reference  
    wing.areas.wetted             = wing.areas.reference  * 2 
    wing.areas.exposed            = wing.areas.reference  * 2  
        
    # add to vehicle 
    vehicle.append_component(wing)   
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Horizontal Tail
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                          = RCAIDE.Components.Wings.Wing()
    wing.tag                      = 'horizontal_tail'  
    wing.aspect_ratio             = 3.04444
    wing.sweeps.quarter_chord     = 17. * Units.degrees
    wing.thickness_to_chord       = 0.12 
    wing.spans.projected          = 2.71805
    wing.chords.root              = 0.94940
    wing.total_length             = 0.94940
    wing.chords.tip               = 0.62731 
    wing.chords.mean_aerodynamic  = 0.809 
    wing.dihedral                 = 20 *Units.degrees
    wing.taper                    = wing.chords.tip / wing.chords.root 
    wing.areas.reference          = 2.14279
    wing.areas.wetted             = 2.14279   * 2
    wing.areas.exposed            = 2.14279   * 2
    wing.twists.root              = 0.0
    wing.twists.tip               = 0.0
    wing.origin                   = [[  5.374 ,0.0 ,  0.596]]
    wing.aerodynamic_center       = [   5.374, 0.0,   0.596] 
    wing.winglet_fraction         = 0.0 
    wing.symmetric                = True    
    
    # add to vehicle 
    vehicle.append_component(wing)        
     
    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Vertical Stabilizer
    #------------------------------------------------------------------------------------------------------------------------------------  
    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'
    
    
    wing.aspect_ratio             = 1.87749
    wing.sweeps.quarter_chord     = 40 * Units.degrees
    wing.thickness_to_chord       = 0.12
    wing.spans.projected          = 1.78701
    wing.chords.root              = 1.79192
    wing.total_length             = 1.79192
    wing.chords.tip               = 0.44578
    wing.taper                    = wing.chords.tip / wing.chords.root 
    wing.chords.mean_aerodynamic  = 0.809  
    wing.areas.reference          = 1.70089
    wing.areas.wetted             = 1.70089 * 2
    wing.areas.exposed            = 1.70089 * 2
    wing.twists.root              = 0.0 * Units.degrees
    wing.twists.tip               = 0.0 * Units.degrees
    wing.origin                   = [[  4.613,0.0 , 0.596]]
    wing.aerodynamic_center       = [0, 0.0, 0]  
    
    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False
    
    wing.dynamic_pressure_ratio  = 1.0
    
    
    # Wing Segments
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.root_chord_percent            = 1.0
    segment.twist                         = 0. * Units.deg 
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 50.0  * Units.degrees  
    segment.thickness_to_chord            = .12   
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 0.41599/wing.spans.projected
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.14447/wing.chords.root   
    segment.dihedral_outboard             = 0. * Units.degrees
    segment.sweeps.quarter_chord          = 40.0 * Units.degrees   
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'segment_2'
    segment.percent_span_location         = 1.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.44578/wing.chords.root 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        
    
    # add to vehicle
    vehicle.append_component(wing) 
      
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################   Fuselage  ############################################################   
    #------------------------------------------------------------------------------------------------------------------------------------ 
    fuselage                                    = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage' 
    fuselage.seats_abreast                      = 2.  
    fuselage.seat_pitch                         = 3.  
    fuselage.fineness.nose                      = 0.88   
    fuselage.fineness.tail                      = 1.13   
    fuselage.lengths.nose                       = 0.5  
    fuselage.lengths.tail                       = 1.5
    fuselage.lengths.cabin                      = 4.46 
    fuselage.lengths.total                      = 6.46
    fuselage.width                              = 5.85 * Units.feet      # change 
    fuselage.heights.maximum                    = 4.65 * Units.feet      # change 
    fuselage.heights.at_quarter_length          = 3.75 * Units.feet      # change 
    fuselage.heights.at_wing_root_quarter_chord = 4.65 * Units.feet      # change 
    fuselage.heights.at_three_quarters_length   = 4.26 * Units.feet      # change 
    fuselage.areas.wetted                       = 236. * Units.feet**2   # change 
    fuselage.areas.front_projected              = 0.14 * Units.feet**2   # change 
    fuselage.effective_diameter                 = 1.276     # change 
    fuselage.differential_pressure              = 0. 
    
    # Segment  
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment() 
    segment.tag                                 = 'segment_0'    
    segment.percent_x_location                  = 0.0 
    segment.percent_z_location                  = 0.     # change  
    segment.height                              = 0.049 
    segment.width                               = 0.032 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_1'   
    segment.percent_x_location                  = 0.10912/fuselage.lengths.total 
    segment.percent_z_location                  = 0.00849
    segment.height                              = 0.481 
    segment.width                               = 0.553 
    fuselage.Segments.append(segment)           
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_2'   
    segment.percent_x_location                  = 0.47804/fuselage.lengths.total
    segment.percent_z_location                  = 0.02874
    segment.height                              = 1.00
    segment.width                               = 0.912 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                            
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  = 0.161  
    segment.percent_z_location                  = 0.04348  
    segment.height                              = 1.41
    segment.width                               = 1.174  
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 0.284 
    segment.percent_z_location                  = 0.05435 
    segment.height                              = 1.62
    segment.width                               = 1.276  
    fuselage.Segments.append(segment)              
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 3.43026/fuselage.lengths.total
    segment.percent_z_location                  = 0.31483/fuselage.lengths.total 
    segment.height                              = 1.409
    segment.width                               = 1.121 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 4.20546/fuselage.lengths.total
    segment.percent_z_location                  = 0.32216/fuselage.lengths.total
    segment.height                              = 1.11
    segment.width                               = 0.833
    fuselage.Segments.append(segment)                  
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 4.99358/fuselage.lengths.total
    segment.percent_z_location                  = 0.37815/fuselage.lengths.total
    segment.height                              = 0.78
    segment.width                               = 0.512 
    fuselage.Segments.append(segment)                  
                                                
    # Segment                                             
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_8'   
    segment.percent_x_location                  = 1.
    segment.percent_z_location                  = 0.55/fuselage.lengths.total
    segment.height                              = 0.195  
    segment.width                               = 0.130 
    fuselage.Segments.append(segment)                   
                                                
    vehicle.append_component(fuselage) 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Booms  ################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------          
    boom                                    = RCAIDE.Components.Fuselages.Fuselage()
    boom.tag                                = 'boom_1r'
    boom.configuration                      = 'boom'  
    boom.origin                             = [[   0.036, 1.950,  1]]  
    boom.seats_abreast                      = 0.  
    boom.seat_pitch                         = 0.0 
    boom.fineness.nose                      = 0.950   
    boom.fineness.tail                      = 1.029   
    boom.lengths.nose                       = 0.2 
    boom.lengths.tail                       = 0.2
    boom.lengths.cabin                      = 4.15
    boom.lengths.total                      = 4.2
    boom.width                              = 0.15 
    boom.heights.maximum                    = 0.15  
    boom.heights.at_quarter_length          = 0.15  
    boom.heights.at_three_quarters_length   = 0.15 
    boom.heights.at_wing_root_quarter_chord = 0.15 
    boom.areas.wetted                       = 0.018
    boom.areas.front_projected              = 0.018 
    boom.effective_diameter                 = 0.15  
    boom.differential_pressure              = 0.  
    boom.symmetric                          = True 
    boom.index                              = 1
    
    # Segment  
    segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment() 
    segment.tag                       = 'segment_1'   
    segment.percent_x_location        = 0.
    segment.percent_z_location        = 0.0 
    segment.height                    = 0.05  
    segment.width                     = 0.05   
    boom.Segments.append(segment)           
    
    # Segment                                   
    segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                       = 'segment_2'   
    segment.percent_x_location        = 0.03
    segment.percent_z_location        = 0. 
    segment.height                    = 0.15 
    segment.width                     = 0.15 
    boom.Segments.append(segment) 
    
    # Segment                                   
    segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                       = 'segment_3'    
    segment.percent_x_location        = 0.97
    segment.percent_z_location        = 0. 
    segment.height                    = 0.15
    segment.width                     = 0.15
    boom.Segments.append(segment)           
    
    # Segment                                  
    segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                       = 'segment_4'   
    segment.percent_x_location        = 1.   
    segment.percent_z_location        = 0.   
    segment.height                    = 0.05   
    segment.width                     = 0.05   
    boom.Segments.append(segment)           
    
    # add to vehicle
    vehicle.append_component(boom)   
    
    # add left long boom 
    boom              = deepcopy(vehicle.fuselages.boom_1r)
    boom.origin[0][1] = -boom.origin[0][1]
    boom.tag          = 'boom_1l' 
    vehicle.append_component(boom)         
     
    # add left long boom 
    boom              = deepcopy(vehicle.fuselages.boom_1r)
    boom.origin       = [[     0.110,    4.891,   1.050]] 
    boom.tag          = 'boom_2r' 
    boom.lengths.total                      = 4.16
    vehicle.append_component(boom)  
     
    # add inner left boom 
    boom              = deepcopy(vehicle.fuselages.boom_1r)
    boom.origin       = [[     0.110, -  4.891,    1.050 ]]   
    boom.lengths.total                      = 4.16
    boom.tag          = 'boom_2l' 
    vehicle.append_component(boom)      
    
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################   Nacelles  ############################################################   
    #------------------------------------------------------------------------------------------------------------------------------------  
    nacelle                = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.tag            = 'rotor_nacelle'
    nacelle.length         = 0.45
    nacelle.diameter       = 0.3
    nacelle.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
    nacelle.flow_through   = False  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.05 
    nac_segment.height             = 0.1
    nac_segment.width              = 0.1
    nacelle.append_segment(nac_segment)    
    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.15 
    nac_segment.height             = 0.2
    nac_segment.width              = 0.2
    nacelle.append_segment(nac_segment)        
    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.25  
    nac_segment.height             = 0.25
    nac_segment.width              = 0.25
    nacelle.append_segment(nac_segment)     
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.25  
    nac_segment.height             = 0.25
    nac_segment.width              = 0.25
    nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.5 
    nac_segment.height             = 0.3
    nac_segment.width              = 0.3
    nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 0.75
    nac_segment.height             = 0.25
    nac_segment.width              = 0.25
    nacelle.append_segment(nac_segment)        
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_8'
    nac_segment.percent_x_location = 1.0
    nac_segment.height             = 0.2
    nac_segment.width              = 0.2
    nacelle.append_segment(nac_segment)      
    
    lift_rotor_nacelle_origins   = [[  -0.073,  1.950, 0.850] ,[ -0.073, -1.950, 0.850],
                           [   4.413,   1.950 ,0.850] ,[   4.413, -1.950, 0.850],
                           [   0.219 ,   4.891 , 0.950] ,[   0.219 , -  4.891 ,0.950],
                           [  4.196 ,   4.891 ,0.950] ,[   4.196, -  4.891 ,0.950]]
    
    for ii in range(8):
        rotor_nacelle          = deepcopy(nacelle)
        rotor_nacelle.tag      = 'rotor_nacelle_' + str(ii+1) 
        rotor_nacelle.origin   = [lift_rotor_nacelle_origins[ii]]
        vehicle.append_component(rotor_nacelle) 
    
    propeller_nacelle                = RCAIDE.Components.Nacelles.Nacelle()
    propeller_nacelle.tag            = 'propeller_nacelle'
    propeller_nacelle.length         = 1.24
    propeller_nacelle.diameter       = 0.4
    propeller_nacelle.orientation_euler_angles  = [0.,0.,0.]    
    propeller_nacelle.flow_through   = False  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.10/propeller_nacelle.length
    nac_segment.height             = 0.2
    nac_segment.width              = 0.2
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.15 /propeller_nacelle.length 
    nac_segment.height             = 0.25
    nac_segment.width              = 0.25
    propeller_nacelle.append_segment(nac_segment)    
    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.2/propeller_nacelle.length  
    nac_segment.height             = 0.3
    nac_segment.width              = 0.3
    propeller_nacelle.append_segment(nac_segment)    
    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.25/propeller_nacelle.length  
    nac_segment.height             = 0.35
    nac_segment.width              = 0.35
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.5/propeller_nacelle.length 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.4
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 0.75/propeller_nacelle.length
    nac_segment.height             = 0.35
    nac_segment.width              = 0.35
    propeller_nacelle.append_segment(nac_segment)        
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_8'
    nac_segment.percent_x_location = 0.98/propeller_nacelle.length
    nac_segment.height             = 0.3
    nac_segment.width              = 0.3
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_9'
    nac_segment.percent_x_location = 1.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    propeller_nacelle.append_segment(nac_segment)     
    
    propeller_nacelle_origins   = [[ 5.583,  1.300 ,    1.092] ,[5.583, - 1.300,     1.092]]
    
    for ii in range(2):
        prop_nacelle          = deepcopy(propeller_nacelle)
        prop_nacelle.tag      = 'propeller_nacelle_' + str(ii+1) 
        prop_nacelle.origin   = [propeller_nacelle_origins[ii]]
        vehicle.append_component(prop_nacelle)    
        
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------     
    sys                            = RCAIDE.Components.Systems.System()
    sys.mass_properties.mass       = 5 # kg   
    vehicle.append_component(sys)   
    
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################  Energy Network  ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    net                              = All_Electric()   
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                              = RCAIDE.Energy.Distributors.Bus_Power_Control_Unit()
    bus.fixed_voltage                = False 
    bus.active_propulsor_groups      = ['forward_propulsor' ,'lift_propulsor']
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Battery   
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat                                                    = RCAIDE.Energy.Storages.Batteries.Lithium_Ion_NMC() 
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 100
    initialize_from_circuit_configuration(bat)  
    bat.module_config.number_of_modules                    = 14 
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module_config.voltage                              = bat.pack.maximum_voltage/bat.module_config.number_of_modules 
    bat.module.geometrtic_configuration.normal_count       = 25
    bat.module.geometrtic_configuration.parallel_count     = 40
    bus.voltage                      =  bat.pack.maximum_voltage  
    bus.batteries.append(bat)                                
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Forward Cruise Propulsor System - Electronic Speed Controller    
    #------------------------------------------------------------------------------------------------------------------------------------  
    propeller_esc                   = RCAIDE.Energy.Distributors.Electronic_Speed_Controller() 
    propeller_esc.efficiency        = 0.95  
    propeller_esc.tag               = 'propeller_esc'  
    propeller_esc.propulsor_group   = 'forward_propulsor'  
    bus.electronic_speed_controllers.append(propeller_esc)    
    
    propeller_esc_2                 = deepcopy(propeller_esc)
    propeller_esc_2.tag             = 'propeller_esc_2'
    propeller_esc_2.propulsor_group = 'forward_propulsor'  
    bus.electronic_speed_controllers.append(propeller_esc_2)   
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Forward Cruise Propulsor System - Propeller    
    #------------------------------------------------------------------------------------------------------------------------------------           
    g                                           = 9.81                                   # gravitational acceleration 
    speed_of_sound                              = 340                                    # speed of sound 
    Drag                                        = estimate_cruise_drag(vehicle,altitude = 1500. * Units.ft,speed= 130.* Units['mph'] ,lift_coefficient = 0.5 ,profile_drag = 0.06)
    Hover_Load                                  = vehicle.mass_properties.takeoff*g      # hover load   
    
    propeller                                   = RCAIDE.Energy.Converters.Propeller()
    propeller.number_of_blades                  = 3
    propeller.tag                               = 'propeller_1'
    propeller.propulsor_group                   = 'forward_propulsor'  
    propeller.tip_radius                        = 1.15  
    propeller.hub_radius                        = 0.1 * propeller.tip_radius  
    propeller.cruise.design_freestream_velocity = 130.* Units['mph'] 
    propeller.cruise.design_tip_mach            = 0.65
    propeller.cruise.design_angular_velocity    = propeller.cruise.design_tip_mach *speed_of_sound/propeller.tip_radius
    propeller.cruise.design_Cl                  = 0.7
    propeller.cruise.design_altitude            = 1500 * Units.feet
    propeller.cruise.design_thrust              = Drag*3/2  
    propeller.rotation                          = 1
    propeller.variable_pitch                    = True  
    airfoil                                     = RCAIDE.Components.Airfoils.Airfoil()
    airfoil.coordinate_file                     = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                         = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                  rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    propeller.append_airfoil(airfoil)          
    propeller.airfoil_polar_stations            = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    propeller                                   = design_propeller(propeller)
                 
    propeller_origins                           = [[  6.583,propeller_nacelle_origins[0][1] , propeller_nacelle_origins[1][2]] ,
                                                   [  6.583,propeller_nacelle_origins[1][1] ,propeller_nacelle_origins[1][2]]]
    propeller.origin                            = [propeller_origins[0]]
    bus.rotors.append(propeller)  
    
    propeller_2          = deepcopy(propeller)
    propeller_2.tag      = 'propeller_2'
    propeller_2.rotation = -1
    propeller_2.origin   = [propeller_origins[1]]
    bus.rotors.append(propeller_2)   
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Forward Cruise Propulsor System - Propeller Motors 
    #------------------------------------------------------------------------------------------------------------------------------------     
    propeller_motor                          = RCAIDE.Energy.Converters.Motor()
    propeller_motor.efficiency               = 0.95
    propeller_motor.tag                      = 'propeller_motor_1'
    propeller_motor.origin                   = propeller.origin
    propeller_motor.nominal_voltage          = bat.pack.maximum_voltage 
    propeller_motor.origin                   = propeller.origin
    propeller_motor.propeller_radius         = propeller.tip_radius
    propeller_motor.propulsor_group          = 'forward_propulsor'
    propeller_motor.no_load_current          = 0.001
    propeller_motor.wing_mounted             = True 
    propeller_motor.wing_tag                 = 'horizontal_tail'
    propeller_motor.rotor_radius             = propeller.tip_radius
    propeller_motor.design_torque            = propeller.cruise.design_torque
    propeller_motor.angular_velocity         = propeller.cruise.design_angular_velocity/propeller_motor.gear_ratio 
    propeller_motor                          = size_optimal_motor(propeller_motor)
    propeller_motor.mass_properties.mass     = nasa_motor(propeller_motor.design_torque) 
    bus.motors.append(propeller_motor) 
    
    propeller_motor_2          = deepcopy(propeller_motor)
    propeller_motor_2.tag      = 'propeller_motor_2' 
    propeller_motor_2.origin   = propeller_2.origin
    bus.motors.append(propeller_motor_2)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Lift Propulsor System - Electronic Speed Controller
    #------------------------------------------------------------------------------------------------------------------------------------    
    lift_rotor_esc                   = RCAIDE.Energy.Distributors.Electronic_Speed_Controller()
    lift_rotor_esc.propulsor_group   = 'lift_propulsor'
    lift_rotor_esc.efficiency        = 0.95
    for i in range(8):
        lift_rotor_ESC          = deepcopy(lift_rotor_esc)
        lift_rotor_ESC.tag      = 'lift_rotor_esc' + str(i + 1)  
        bus.electronic_speed_controllers.append(lift_rotor_ESC) 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Lift Propulsor System - Lift Rotors   
    #------------------------------------------------------------------------------------------------------------------------------------        
    lift_rotor                                   = RCAIDE.Energy.Converters.Lift_Rotor()  
    lift_rotor.tip_radius                        = 2.8/2
    lift_rotor.hub_radius                        = 0.15 * lift_rotor.tip_radius   
    lift_rotor.number_of_blades                  = 3     
    lift_rotor.hover.design_altitude             = 40 * Units.feet  
    lift_rotor.hover.design_thrust               = Hover_Load/8
    lift_rotor.propulsor_group                   = 'lift_propulsor'
    lift_rotor.hover.design_freestream_velocity  = np.sqrt(lift_rotor.hover.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    lift_rotor.oei.design_altitude               = 40 * Units.feet  
    lift_rotor.oei.design_thrust                 = Hover_Load/6  
    lift_rotor.oei.design_freestream_velocity    = np.sqrt(lift_rotor.oei.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    airfoil                                      = RCAIDE.Components.Airfoils.Airfoil()   
    airfoil.coordinate_file                      = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                          =[rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                   rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                    rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    lift_rotor.append_airfoil(airfoil)               
    lift_rotor.airfoil_polar_stations           = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    lift_rotor                                  = design_lift_rotor(lift_rotor)   
    
    lift_rotor_origins                     = [[ -0.073, 1.950,  1.2] ,  [  -0.073, -1.950 ,  1.2],
                                              [ 4.440 , 1.950 ,  1.2] ,[ 4.440 , -1.950,  1.2],
                                              [ 0.219 ,  4.891 , 1.2]  ,[ 0.219 , - 4.891 , 1.2],
                                              [ 4.196 ,  4.891 , 1.2] ,[4.196, - 4.891 , 1.2]]
    
    # Appending rotors with different origins
    rotations = [-1,1,-1,1,-1,1,-1,1] 
    angle_offsets        = np.random.rand(8)*(np.pi)    
    for ii in range(8):
        lift_rotor                        = deepcopy(lift_rotor)
        lift_rotor.tag                    = 'lift_rotor_' + str(ii+1)
        lift_rotor.rotation               = rotations[ii]
        lift_rotor.origin                 = [lift_rotor_origins[ii]]
        lift_rotor.phase_offset_angle     = angle_offsets[ii]
        bus.rotors.append(lift_rotor)    
        
    
    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Propulsor System - Lift Rotor Motors 
    #------------------------------------------------------------------------------------------------------------------------------------  
    lift_rotor_motor                         = RCAIDE.Energy.Converters.Motor()
    lift_rotor_motor.efficiency              = 0.9
    lift_rotor_motor.nominal_voltage         = bat.pack.maximum_voltage*3/4  
    lift_rotor_motor.origin                  = lift_rotor.origin 
    lift_rotor_motor.propeller_radius        = lift_rotor.tip_radius
    lift_rotor_motor.propulsor_group         = 'lift_propulsor'   
    lift_rotor_motor.no_load_current         = 0.01  
    lift_rotor_motor.wing_mounted            = True 
    lift_rotor_motor.wing_tag                = 'main_wing'
    lift_rotor_motor.rotor_radius            = lift_rotor.tip_radius
    lift_rotor_motor.design_torque           = lift_rotor.hover.design_torque
    lift_rotor_motor.angular_velocity        = lift_rotor.hover.design_angular_velocity/lift_rotor_motor.gear_ratio  
    lift_rotor_motor                         = size_optimal_motor(lift_rotor_motor)
    lift_rotor_motor.mass_properties.mass    = nasa_motor(lift_rotor_motor.design_torque)    
    
    # Appending motors with different origins
    for i in range(8):
        lr_motor           = deepcopy(lift_rotor_motor)
        lr_motor.tag       = 'lift_rotor_motor_' + str(i+1)
        lift_rotor.origin  = [lift_rotor_origins[ii]]
        bus.motors.append(lr_motor) 
    
    # append bus   
    net.busses.append(bus)
    
    #------------------------------------------------------------------------------------------------------------------------------------           
    # Payload 
    #------------------------------------------------------------------------------------------------------------------------------------  
    payload                        = RCAIDE.Energy.Peripherals.Avionics()
    payload.power_draw             = 10. # Watts 
    payload.mass_properties.mass   = 1.0 * Units.kg
    net.payload                    = payload
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------  
    avionics                       = RCAIDE.Energy.Peripherals.Avionics()
    avionics.power_draw            = 20. # Watts  
    avionics.mass_properties.mass  = 1.0 * Units.kg
    net.avionics                   = avionics  
     
    # append energy network 
    vehicle.append_energy_network(net)
    
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------      
    settings = Data()
    converge_evtol_weight(vehicle,settings,contingency_factor = 1.0) 
    breakdown = compute_weight(vehicle,settings,contingency_factor = 1.0 )
    print(breakdown)
    
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
    base_config.networks.all_electric.busses.bus.active_propulsor_groups = ['forward_propulsor' ,'lift_propulsor']    
    configs.append(base_config)


    forward_config = RCAIDE.Components.Configs.Config(vehicle)
    forward_config.tag = 'forward_flight' 
    forward_config.networks.all_electric.busses.bus.active_propulsor_groups = ['forward_propulsor'] 
    configs.append(forward_config) 


    transition_config = RCAIDE.Components.Configs.Config(vehicle)
    transition_config.tag = 'transition_flight' 
    transition_config.networks.all_electric.busses.bus.active_propulsor_groups = ['forward_propulsor' ,'lift_propulsor']    
    configs.append(transition_config)
    

    vertical_config = RCAIDE.Components.Configs.Config(vehicle)
    vertical_config.tag = 'vertical_flight'  
    vertical_config.networks.all_electric.busses.bus.active_propulsor_groups = ['lift_propulsor'] 
    configs.append(vertical_config)  


    descent_config = RCAIDE.Components.Configs.Config(vehicle)
    descent_config.tag = 'descent' 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_1.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_2.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_3.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_4.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_5.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_6.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_7.pitch_command           =-5 * Units.degrees 
    base_config.networks.all_electric.busses.bus.rotors.lift_rotor_8.pitch_command           =-5 * Units.degrees  
    descent_config.networks.all_electric.busses.bus.active_propulsor_groups = ['forward_propulsor'] 
    configs.append(descent_config)  
    
    
    # done!
    return configs

  
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
    weights         = RCAIDE.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics          = RCAIDE.Analyses.Aerodynamics.Fidelity_Zero() 
    aerodynamics.geometry = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)   

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Analyses.Energy.Energy()
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

    # done!
    return analyses    

# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def mission_setup(analyses): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission' 
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()   
    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
     
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                          = 0.0  * Units.ft  
    segment.altitude_end                            = 200.  * Units.ft   
    segment.initial_battery_state_of_charge         = 1.0 
    segment.climb_rate                              = 500. * Units['ft/min']  
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight )   
    segment.altitude                                = 200.  * Units.ft           
    segment.air_speed_start                         = 500. * Units['ft/min']
    segment.air_speed_end                           = 0.75 * Vstall
    segment.acceleration                            = 1.5
    segment.pitch_initial                           = 0.0 * Units.degrees
    segment.pitch_final                             = 2.  * Units.degrees    
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                     = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)    
    segment.altitude_start                          = 200.0 * Units.ft   
    segment.altitude_end                            = 500.0 * Units.ft   
    segment.air_speed_start                         = 0.75   * Vstall
    segment.climb_angle                             = 3     * Units.degrees   
    segment.acceleration                            = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                           = 2.    * Units.degrees 
    segment.pitch_final                             = 7.    * Units.degrees    
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment) 
  
    # ------------------------------------------------------------------
    #   First Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_1"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft   
    segment.altitude_end                             = 1000. * Units.ft   
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_start                          = 104.1  * Units['mph']  
    segment.air_speed_end                            = 105.  * Units['mph']        
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #  Second Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_2"  
    segment.analyses.extend( analyses.forward_flight) 
    segment.altitude_start                           = 1000.0 * Units.ft   
    segment.altitude_end                             = 1500. * Units.ft   
    segment.climb_rate                               = 300.  * Units['ft/min']
    segment.air_speed_start                          = 105.  * Units['mph']  
    segment.air_speed_end                            = 110.  * Units['mph']        
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])     
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Third Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise"  
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 1500.0 * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 50 *Units.nmi 
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])      
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                           = 1500.0 * Units.ft  
    segment.altitude_end                             = 1000. * Units.ft  
    segment.climb_rate                               = -500.  * Units['ft/min']
    segment.air_speed_start                          = 110.  * Units['mph']  
    segment.air_speed_end                            = Vstall     
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])      
    mission.append_segment(segment)  
     
            
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                      = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude_start                           = 1000.0 * Units.ft   
    segment.altitude_end                             = 300.0 * Units.ft   
    segment.climb_angle                              = 3 * Units.degrees
    segment.acceleration                             = -0.25 * Units['m/s/s']    
    segment.pitch_initial                            = 4.3  * Units.degrees  
    segment.pitch_final                              = 7. * Units.degrees            
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 300.  * Units.ft     
    segment.air_speed_start                         = 36.5 * Units['mph'] 
    segment.air_speed_end                           = 300. * Units['ft/min'] 
    segment.acceleration                            = -0.25 * Units['m/s/s']    
    segment.pitch_initial                           = 7.  * Units.degrees  
    segment.pitch_final                             = 7. * Units.degrees      
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Hover.Descent(base_segment)
    segment.tag                                     = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                          = 300.0 * Units.ft   
    segment.altitude_end                            = 0.   * Units.ft  
    segment.descent_rate                            = 300. * Units['ft/min']      
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)  
  

    # ------------------------------------------------------------------
    #  Charge Segment: 
    # ------------------------------------------------------------------  
    # Charge Model 
    segment                                         = Segments.Ground.Battery_Recharge(base_segment)     
    segment.tag                                     = 'Recharge' 
    segment.time                                    = 1 * Units.hr
    segment.analyses.extend(analyses.base)              
    segment = analyses.base.energy.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)    
    return mission 


def missions_setup(mission): 
 
    missions         = RCAIDE.Analyses.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  


# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------  
def plot_results(results,run_noise_model,save_figure_flag):  
    
    # Plots fligh conditions 
    plot_flight_conditions(results) 
    
    # Plot arcraft trajectory
    plot_flight_trajectory(results)
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
     
    # Plot Aircraft Stability
    plot_stability_coefficients(results) 
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results) 
    plot_battery_health_conditions(results)
    plot_battery_cell_conditions(results) 
    plot_battery_degradation(results) 
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
    plot_disc_and_power_loading(results)
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_motor_and_rotor_efficiencies(results) 
    
    # Plot Battery Degradation  
    plot_battery_degradation(results)   
     
    return
 
if __name__ == '__main__': 
    main()    
    plt.show() 