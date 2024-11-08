'''

The script below documents how to set up and plot the results of polar analysis of full aircraft configuration 

''' 

# ----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------- 
import RCAIDE
from RCAIDE.Framework.Core import Units , Data   
from RCAIDE.Library.Methods.Propulsors.Turbojet_Propulsor          import design_turbojet
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor          import design_turbofan  
from RCAIDE.Library.Methods.Stability.Center_of_Gravity            import compute_component_centers_of_gravity 
from RCAIDE.Library.Methods.Geometry.Planform                      import segment_properties, wing_segmented_planform 
from RCAIDE.Library.Plots                                          import *        

import matplotlib.cm as cm   
import numpy as np
from copy import deepcopy
import matplotlib.pyplot  as plt
import os

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main(): 

    vehicle                            = Boeing_737_vehicle_setup()  
    Mach_range                         = np.linspace(0.1, 0.9, 20)
    fig_title                          =  'Subsonic'
    subsonic_aero_coefficient_test(vehicle,Mach_range,fig_title )

    vehicle                            = Concorde_vehicle_setup()  
    Mach_range                         = np.linspace(0.2,3,20)
    fig_title                          = 'Supersonic'
    supersonic_aero_coefficient_test(vehicle,Mach_range,fig_title )    

    control_surface_deflection_angles  = np.linspace(0,30,7)*Units.degrees     
    subsonic_control_surface_test(control_surface_deflection_angles)
    
    return 
    
def subsonic_aero_coefficient_test(vehicle,Mach_range,fig_title ):
    altitude                           = 5000.0*Units.feet
    alpha_range                        = np.linspace(-2,10,13)*Units.degrees 
    MAC                                = vehicle.wings['main_wing'].chords.mean_aerodynamic 
      

    #------------------------------------------------------------------------
    # setup flight conditions
    #------------------------------------------------------------------------   
    AoA_range  = np.atleast_2d(alpha_range).T 
    MAC            = vehicle.wings['main_wing'].chords.mean_aerodynamic 
    atmosphere     = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data      = atmosphere.compute_values(altitude=altitude)
    P              = atmo_data.pressure[0]
    T              = atmo_data.temperature[0]
    rho            = atmo_data.density[0]  
    a              = atmo_data.speed_of_sound[0]
    mu             = atmo_data.dynamic_viscosity[0]  
       
    # -----------------------------------------------------------------
    # Evaluate Without Surrogate
    # ----------------------------------------------------------------- 
    
    state                                         = RCAIDE.Framework.Mission.Common.State()
    state.conditions                              = RCAIDE.Framework.Mission.Common.Results() 
    state.conditions.freestream.density           = rho * np.ones_like(AoA_range)
    state.conditions.freestream.dynamic_viscosity = mu  * np.ones_like(AoA_range)
    state.conditions.freestream.temperature       = T   * np.ones_like(AoA_range)
    state.conditions.freestream.pressure          = P   * np.ones_like(AoA_range)
    state.conditions.aerodynamics.angles.alpha    = AoA_range 
    state.conditions.aerodynamics.angles.beta     = AoA_range *0  
    state.conditions.freestream.u                 = AoA_range *0       
    state.conditions.freestream.v                 = AoA_range *0       
    state.conditions.freestream.w                 = AoA_range *0       
    state.conditions.static_stability.roll_rate   = AoA_range *0       
    state.conditions.static_stability.pitch_rate  = AoA_range *0 
    state.conditions.static_stability.yaw_rate    = AoA_range *0  

    CL_no_surrogate = np.zeros((len(AoA_range),len(Mach_range)))
    CD_no_surrogate = np.zeros((len(AoA_range),len(Mach_range)))
    CL_surrogate    = np.zeros((len(AoA_range),len(Mach_range)))  
    CD_surrogate    = np.zeros((len(AoA_range),len(Mach_range))) 

    state.analyses                                         =  Data()
    no_sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    no_sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    no_sur_aerodynamics.settings.use_surrogate             = False 
    no_sur_aerodynamics.vehicle                           = vehicle
    no_sur_aerodynamics.settings.model_fuselage            = True   
    no_sur_aerodynamics.initialize()
    state.analyses.aerodynamics = no_sur_aerodynamics           

    # ---------------------------------------------------------------------------------------
    # Evaluate With Surrogate
    # ---------------------------------------------------------------------------------------

    state.analyses                                      =  Data()
    sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    sur_aerodynamics.settings.use_surrogate             = True 
    sur_aerodynamics.vehicle                           = vehicle
    sur_aerodynamics.settings.model_fuselage            = True   
    sur_aerodynamics.initialize() 
    state.analyses.aerodynamics = sur_aerodynamics 
    
    
    for i in range (len(Mach_range)):  
        state.conditions.freestream.mach_number       = Mach_range[i] * np.ones_like(AoA_range)
        state.conditions.freestream.velocity          = Mach_range[i] * a   * np.ones_like(AoA_range) 
        state.conditions.freestream.reynolds_number   =  (Mach_range[i] * a*rho*MAC)/mu   * np.ones_like(AoA_range)
        state.conditions.aerodynamics.angles.alpha    = AoA_range 
     
        _                      = no_sur_aerodynamics.evaluate(state)
        CL_no_surrogate[:,i]   = state.conditions.aerodynamics.coefficients.lift.total[:, 0]
        CD_no_surrogate[:,i]   = state.conditions.aerodynamics.coefficients.drag.total[:, 0]
    
        # ---------------------------------------------------------------------------------------
        # Evaluate With Surrogate
        # --------------------------------------------------------------------------------------- 
        _                      = sur_aerodynamics.evaluate(state)        
        CL_surrogate[:,i]      = state.conditions.aerodynamics.coefficients.lift.total[:, 0]
        CD_surrogate[:,i]      = state.conditions.aerodynamics.coefficients.drag.total[:, 0] 
         

    #------------------------------------------------------------------------
    # setup figures
    #------------------------------------------------------------------------
    fig_1 = plt.figure(fig_title + " Lift Coefficients") 
    fig_2 = plt.figure(fig_title + " Drag Coefficients") 
    fig_1.set_size_inches(12,6)
    fig_2.set_size_inches(12,6)
    axis_1 = fig_1.add_subplot(1, 2, 1, projection='3d')
    axis_2 = fig_1.add_subplot(1, 2, 2, projection='3d')
    axis_3 = fig_2.add_subplot(1, 2, 1, projection='3d')
    axis_4 = fig_2.add_subplot(1, 2, 2, projection='3d')
    # Make data.
    
    X, Y = np.meshgrid(Mach_range, AoA_range)
    axis_1.set_title('$C_L$ Surrogate') 
    axis_2.set_title('$C_L$ No Surrogate')    
    axis_3.set_title('$C_D$ Surrogate') 
    axis_4.set_title('$C_D$ No Surrogate')    
    surf = axis_1.plot_surface(X, Y/Units.degree, CL_surrogate   , cmap=cm.jet,linewidth=0, antialiased=False)
    surf = axis_2.plot_surface(X, Y/Units.degree, CL_no_surrogate, cmap=cm.jet,linewidth=0, antialiased=False)   
    surf = axis_3.plot_surface(X, Y/Units.degree, CD_surrogate   , cmap=cm.jet,linewidth=0, antialiased=False)
    surf = axis_4.plot_surface(X, Y/Units.degree, CD_no_surrogate, cmap=cm.jet,linewidth=0, antialiased=False)    
       
    axis_1.set_ylabel('AoA') 
    axis_2.set_ylabel('AoA') 
    axis_3.set_ylabel('AoA') 
    axis_4.set_ylabel('AoA') 
    axis_1.set_xlabel('Mach') 
    axis_2.set_xlabel('Mach')  
    axis_3.set_xlabel('Mach') 
    axis_4.set_xlabel('Mach')
    
    plt.tight_layout()
    return



def supersonic_aero_coefficient_test(vehicle,Mach_range,fig_title ):
    altitude                           = 5000.0*Units.feet
    alpha_range                        = np.linspace(-2,10,13)*Units.degrees 
    MAC                                = vehicle.wings['main_wing'].chords.mean_aerodynamic 
      

    #------------------------------------------------------------------------
    # setup flight conditions
    #------------------------------------------------------------------------   
    AoA_range  = np.atleast_2d(alpha_range).T  
    atmosphere     = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data      = atmosphere.compute_values(altitude=altitude)
    P              = atmo_data.pressure[0]
    T              = atmo_data.temperature[0]
    rho            = atmo_data.density[0]  
    a              = atmo_data.speed_of_sound[0]
    mu             = atmo_data.dynamic_viscosity[0]  
       
    # -----------------------------------------------------------------
    # Evaluate Without Surrogate
    # ----------------------------------------------------------------- 
    
    state                                         = RCAIDE.Framework.Mission.Common.State()
    state.conditions                              = RCAIDE.Framework.Mission.Common.Results() 
    state.conditions.freestream.density           = rho * np.ones_like(AoA_range)
    state.conditions.freestream.dynamic_viscosity = mu  * np.ones_like(AoA_range)
    state.conditions.freestream.temperature       = T   * np.ones_like(AoA_range)
    state.conditions.freestream.pressure          = P   * np.ones_like(AoA_range)
    state.conditions.aerodynamics.angles.alpha    = AoA_range 
    state.conditions.aerodynamics.angles.beta     = AoA_range *0  
    state.conditions.freestream.u                 = AoA_range *0       
    state.conditions.freestream.v                 = AoA_range *0       
    state.conditions.freestream.w                 = AoA_range *0       
    state.conditions.static_stability.roll_rate   = AoA_range *0       
    state.conditions.static_stability.pitch_rate  = AoA_range *0 
    state.conditions.static_stability.yaw_rate    = AoA_range *0  

    CL_no_surrogate = np.zeros((len(AoA_range),len(Mach_range)))
    CD_no_surrogate = np.zeros((len(AoA_range),len(Mach_range)))
    CL_surrogate    = np.zeros((len(AoA_range),len(Mach_range)))  
    CD_surrogate    = np.zeros((len(AoA_range),len(Mach_range)))
    
 
    state.analyses                                         =  Data()
    no_sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    no_sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    no_sur_aerodynamics.settings.use_surrogate             = False 
    no_sur_aerodynamics.vehicle                           = vehicle
    no_sur_aerodynamics.settings.model_fuselage            = True   
    no_sur_aerodynamics.initialize()
    state.analyses.aerodynamics = no_sur_aerodynamics           
    # ---------------------------------------------------------------------------------------
    # Evaluate With Surrogate
    # ---------------------------------------------------------------------------------------

    state.analyses                                      =  Data()
    sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    sur_aerodynamics.settings.use_surrogate             = True 
    sur_aerodynamics.vehicle                           = vehicle
    sur_aerodynamics.settings.model_fuselage            = True   
    sur_aerodynamics.initialize() 
    state.analyses.aerodynamics = sur_aerodynamics 
    
    
    for i in range (len(Mach_range)):  
        state.conditions.freestream.mach_number       = Mach_range[i] * np.ones_like(AoA_range)
        state.conditions.freestream.velocity          = Mach_range[i] * a   * np.ones_like(AoA_range) 
        state.conditions.freestream.reynolds_number   =  (Mach_range[i] * a*rho*MAC)/mu   * np.ones_like(AoA_range)
        state.conditions.aerodynamics.angles.alpha    = AoA_range 
     
        _                      = no_sur_aerodynamics.evaluate(state)
        CL_no_surrogate[:,i]   = state.conditions.aerodynamics.coefficients.lift.total[:, 0]
        CD_no_surrogate[:,i]   = state.conditions.aerodynamics.coefficients.drag.total[:, 0]
    
        # ---------------------------------------------------------------------------------------
        # Evaluate With Surrogate
        # --------------------------------------------------------------------------------------- 
        state.analyses.aerodynamics = sur_aerodynamics
    
        _                      = sur_aerodynamics.evaluate(state)        
        CL_surrogate[:,i]      = state.conditions.aerodynamics.coefficients.lift.total[:, 0]
        CD_surrogate[:,i]      = state.conditions.aerodynamics.coefficients.drag.total[:, 0] 
         

    #------------------------------------------------------------------------
    # setup figures
    #------------------------------------------------------------------------
    fig_1 = plt.figure(fig_title + " Lift Coefficients") 
    fig_2 = plt.figure(fig_title + " Drag Coefficients") 
    fig_1.set_size_inches(12,6)
    fig_2.set_size_inches(12,6)
    axis_1 = fig_1.add_subplot(1, 2, 1, projection='3d')
    axis_2 = fig_1.add_subplot(1, 2, 2, projection='3d')
    axis_3 = fig_2.add_subplot(1, 2, 1, projection='3d')
    axis_4 = fig_2.add_subplot(1, 2, 2, projection='3d')
    # Make data.
    
    X, Y = np.meshgrid(Mach_range, AoA_range)
    axis_1.set_title('$C_L$ Surrogate') 
    axis_2.set_title('$C_L$ No Surrogate')    
    axis_3.set_title('$C_D$ Surrogate') 
    axis_4.set_title('$C_D$ No Surrogate')    
    surf = axis_1.plot_surface(X, Y/Units.degree, CL_surrogate   , cmap=cm.jet,linewidth=0, antialiased=False)
    surf = axis_2.plot_surface(X, Y/Units.degree, CL_no_surrogate, cmap=cm.jet,linewidth=0, antialiased=False)   
    surf = axis_3.plot_surface(X, Y/Units.degree, CD_surrogate   , cmap=cm.jet,linewidth=0, antialiased=False)
    surf = axis_4.plot_surface(X, Y/Units.degree, CD_no_surrogate, cmap=cm.jet,linewidth=0, antialiased=False)    
       
    axis_1.set_ylabel('AoA') 
    axis_2.set_ylabel('AoA') 
    axis_3.set_ylabel('AoA') 
    axis_4.set_ylabel('AoA') 
    axis_1.set_xlabel('Mach') 
    axis_2.set_xlabel('Mach')  
    axis_3.set_xlabel('Mach') 
    axis_4.set_xlabel('Mach')
    
    plt.tight_layout()
    return
    
def subsonic_control_surface_test(control_surface_deflection_angles): 
    vehicle                            = Boeing_737_vehicle_setup() 
    altitude                           = 5000.0*Units.feet
    Mach                               = 0.5
    alpha_range                        = np.linspace(-2,10,13)*Units.degrees 
    MAC                                = vehicle.wings['main_wing'].chords.mean_aerodynamic
    S_ref                              = vehicle.reference_area   
    num_def                            = len(control_surface_deflection_angles)
    
    #------------------------------------------------------------------------
    # setup figures
    #------------------------------------------------------------------------
    fig = plt.figure('Control_Surfaces_Test')
    fig.set_size_inches(12,6)
    fig.suptitle('Boeing 737: $S_{ref}$ = ' + str(round(S_ref,2)) + ' MAC =' + str(round(MAC,2)))
    axes1 = fig.add_subplot(1,3,1)
    axes1.set_title('$C_L$ vs AoA')
    axes1.set_ylabel('Coefficient of Lift')
    axes1.set_xlabel('Angle of Attack (degrees)') 

    axes2 = fig.add_subplot(1,3,2)
    axes2.set_title('$C_{D}$ vs. AoA')
    axes2.set_ylabel('Coefficient of Drag')
    axes2.set_xlabel('Angle of Attack (degrees)') 

    axes3 = fig.add_subplot(1,3,3)
    axes3.set_title('$C_{D}$ vs. $C_{L}$')
    axes3.set_ylabel('Coefficient of Lift')
    axes3.set_xlabel('Coefficient of Drag')  
    
    linecolor_1 = cm.jet(np.linspace(0, 1,num_def))   
    linecolor_2 = cm.jet(np.linspace(0, 1,num_def)) 
    linestyle_1 = '-' 
    marker_1    = None 
    marker_2    = 's'

    #------------------------------------------------------------------------
    # setup flight conditions
    #------------------------------------------------------------------------   
    AoA_range  = np.atleast_2d(alpha_range).T 
    MAC            = vehicle.wings['main_wing'].chords.mean_aerodynamic 
    atmosphere     = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data      = atmosphere.compute_values(altitude=altitude)
    P              = atmo_data.pressure[0]
    T              = atmo_data.temperature[0]
    rho            = atmo_data.density[0]  
    a              = atmo_data.speed_of_sound[0]
    mu             = atmo_data.dynamic_viscosity[0] 
    V              = a*Mach
    re             = (V*rho*MAC)/mu  
       
    # -----------------------------------------------------------------
    # Evaluate Without Surrogate
    # ----------------------------------------------------------------- 
    
    state                                         = RCAIDE.Framework.Mission.Common.State()
    state.conditions                              = RCAIDE.Framework.Mission.Common.Results() 
    state.conditions.freestream.mach_number       = Mach  * np.ones_like(AoA_range)
    state.conditions.freestream.density           = rho * np.ones_like(AoA_range)
    state.conditions.freestream.dynamic_viscosity = mu  * np.ones_like(AoA_range)
    state.conditions.freestream.temperature       = T   * np.ones_like(AoA_range)
    state.conditions.freestream.pressure          = P   * np.ones_like(AoA_range)
    state.conditions.freestream.reynolds_number   = re  * np.ones_like(AoA_range)
    state.conditions.freestream.velocity          = V   * np.ones_like(AoA_range)
    state.conditions.aerodynamics.angles.alpha    = AoA_range 
    state.conditions.aerodynamics.angles.beta     = AoA_range *0  
    state.conditions.freestream.u                 = AoA_range *0       
    state.conditions.freestream.v                 = AoA_range *0       
    state.conditions.freestream.w                 = AoA_range *0       
    state.conditions.static_stability.roll_rate   = AoA_range *0       
    state.conditions.static_stability.pitch_rate  = AoA_range *0 
    state.conditions.static_stability.yaw_rate    = AoA_range *0
    

    state.analyses                                         =  Data()
    no_sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    no_sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    no_sur_aerodynamics.settings.use_surrogate             = False 
    no_sur_aerodynamics.vehicle                           = vehicle
    no_sur_aerodynamics.settings.model_fuselage            = True   
    no_sur_aerodynamics.initialize()
    state.analyses.aerodynamics = no_sur_aerodynamics           

    # ---------------------------------------------------------------------------------------
    # Evaluate With Surrogate
    # ---------------------------------------------------------------------------------------

    state.analyses                                      =  Data()
    sur_aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    sur_aerodynamics.settings.fuselage_lift_correction  = 1. 
    sur_aerodynamics.settings.use_surrogate             = True 
    sur_aerodynamics.vehicle                           = vehicle
    sur_aerodynamics.settings.model_fuselage            = True   
    sur_aerodynamics.initialize() 
    state.analyses.aerodynamics = sur_aerodynamics 
    
     
    for i in range (num_def):  
        state.conditions.aerodynamics.angles.alpha    = AoA_range 
        no_sur_aerodynamics.vehicle.wings.main_wing.control_surfaces.flap.deflection = control_surface_deflection_angles[i]
        
        _                 = no_sur_aerodynamics.evaluate(state)
        CL_no_surrogate   = state.conditions.aerodynamics.coefficients.lift.total
        CD_no_surrogate   = state.conditions.aerodynamics.coefficients.drag.total
    
        # ---------------------------------------------------------------------------------------
        # Evaluate With Surrogate
        # --------------------------------------------------------------------------------------- 

        sur_aerodynamics.vehicle.wings.main_wing.control_surfaces.flap.deflection = control_surface_deflection_angles[i]    
        _                 = sur_aerodynamics.evaluate(state)        
        CL_surrogate      = state.conditions.aerodynamics.coefficients.lift.total
        CD_surrogate      = state.conditions.aerodynamics.coefficients.drag.total 
         
        axes1.plot(AoA_range/Units.degrees,CL_no_surrogate,linestyle = linestyle_1, color = linecolor_1[i], marker = marker_1) 
        axes1.scatter(AoA_range/Units.degrees,CL_surrogate , c = linecolor_2[i], s = 50, marker  = marker_2)  

        line_label =  vehicle.wings.main_wing.tag + 'flap : ' + str(round( control_surface_deflection_angles[i]/Units.degrees,3)) + '$\degree$ defl.'  
        axes2.plot(AoA_range/Units.degrees,CD_no_surrogate,linestyle = linestyle_1, color = linecolor_1[i], marker = marker_1,label = line_label) 
        axes2.scatter(AoA_range/Units.degrees,CD_surrogate , c = linecolor_2[i], s = 50,marker = marker_2) 
         
        axes3.plot(CD_no_surrogate,CL_no_surrogate,linestyle = linestyle_1, color = linecolor_1[i], marker = marker_1) 
        axes3.scatter(CD_surrogate,CL_surrogate,  c= linecolor_2[i],  s = 50,marker = marker_2)
         
    axes2.legend(loc='upper left', prop={'size': 8}) 
    plt.tight_layout()
    return 
 
 
# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def Boeing_737_vehicle_setup():  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Boeing_737-800'    
    vehicle.mass_properties.max_takeoff               = 79015.8 * Units.kilogram  
    vehicle.mass_properties.takeoff                   = 79015.8 * Units.kilogram    
    vehicle.mass_properties.operating_empty           = 62746.4 * Units.kilogram  
    vehicle.mass_properties.max_zero_fuel             = 62732.0 * Units.kilogram 
    vehicle.mass_properties.cargo                     = 10000.  * Units.kilogram  
    vehicle.envelope.ultimate_load                    = 3.75
    vehicle.envelope.limit_load                       = 2.5 
    vehicle.reference_area                            = 124.862 * Units['meters**2']   
    vehicle.passengers                                = 170
    vehicle.systems.control                           = "fully powered" 
    vehicle.systems.accessories                       = "medium range"
 
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################################### Landing Gear ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------
    landing_gear                    = RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                = "main_landing_gear" 
    landing_gear.main_tire_diameter = 1.12000 * Units.m
    landing_gear.nose_tire_diameter = 0.6858 * Units.m
    landing_gear.main_strut_length  = 1.8 * Units.m
    landing_gear.nose_strut_length  = 1.3 * Units.m
    landing_gear.main_units         = 2    # Number of main landing gear
    landing_gear.nose_units         = 1    # Number of nose landing gear
    landing_gear.main_wheels        = 2    # Number of wheels on the main landing gear
    landing_gear.nose_wheels        = 2    # Number of wheels on the nose landing gear      
    vehicle.landing_gear            = landing_gear

    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.aspect_ratio                     = 10.18
    wing.sweeps.quarter_chord             = 25 * Units.deg
    wing.thickness_to_chord               = 0.1
    wing.taper                            = 0.1 
    wing.spans.projected                  = 34.32 
    wing.chords.root                      = 7.760 * Units.meter
    wing.chords.tip                       = 0.782 * Units.meter
    wing.chords.mean_aerodynamic          = 4.235 * Units.meter 
    wing.areas.reference                  = 124.862
    wing.areas.wetted                     = 225.08 
    wing.twists.root                      = 4.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[13.61,0,-0.5]]
    wing.aerodynamic_center               = [0,0,0] 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.dynamic_pressure_ratio           = 1.0


    # Wing Segments
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator + '..' + separator  
    root_airfoil.coordinate_file          = rel_path  + 'Airfoils' + separator + 'B737a.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 4. * Units.deg
    segment.root_chord_percent            = 1.
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 2.5 * Units.degrees
    segment.sweeps.quarter_chord          = 28.225 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    yehudi_airfoil                        = RCAIDE.Library.Components.Airfoils.Airfoil()
    yehudi_airfoil.coordinate_file        = rel_path+ 'Airfoils' + separator + 'B737b.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Yehudi'
    segment.percent_span_location         = 0.324
    segment.twist                         = 0.047193 * Units.deg
    segment.root_chord_percent            = 0.5
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 25. * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(yehudi_airfoil)
    wing.append_segment(segment)

    mid_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    mid_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737c.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Section_2'
    segment.percent_span_location         = 0.963
    segment.twist                         = 0.00258 * Units.deg
    segment.root_chord_percent            = 0.220
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 56.75 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(mid_airfoil)
    wing.append_segment(segment)

    tip_airfoil                           =  RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737d.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees
    segment.root_chord_percent            = 0.10077
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = .1
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)    

    # control surfaces -------------------------------------------
    slat                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Slat()
    slat.tag                      = 'slat'
    slat.span_fraction_start      = 0.2
    slat.span_fraction_end        = 0.963
    slat.deflection               = 0.0 * Units.degrees
    slat.chord_fraction           = 0.075
    wing.append_control_surface(slat)

    flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.7
    flap.deflection               = 0.0 * Units.degrees
    flap.configuration_type       = 'double_slotted'
    flap.chord_fraction           = 0.30
    wing.append_control_surface(flap)

    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.963
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16
    wing.append_control_surface(aileron)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing     = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 4.99
    wing.sweeps.quarter_chord    = 28.2250 * Units.deg  
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.3333  
    wing.spans.projected         = 14.4 
    wing.chords.root             = 4.2731 
    wing.chords.tip              = 1.4243 
    wing.chords.mean_aerodynamic = 8.0 
    wing.areas.reference         = 41.49
    wing.areas.exposed           = 59.354    # Exposed area of the horizontal tail
    wing.areas.wetted            = 71.81     # Wetted area of the horizontal tail
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees 
    wing.origin                  = [[33.02,0,1.466]]
    wing.aerodynamic_center      = [0,0,0] 
    wing.vertical                = False
    wing.symmetric               = True 
    wing.dynamic_pressure_ratio  = 0.9


    # Wing Segments
    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 8.63 * Units.degrees
    segment.sweeps.quarter_chord   = 28.2250  * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 0.3333               
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        

    # control surfaces -------------------------------------------
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.09
    elevator.span_fraction_end     = 0.92
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 1.98865
    wing.sweeps.quarter_chord    = 31.2  * Units.deg   
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.1183

    wing.spans.projected         = 8.33
    wing.total_length            = wing.spans.projected 
    
    wing.chords.root             = 10.1 
    wing.chords.tip              = 1.20 
    wing.chords.mean_aerodynamic = 4.0

    wing.areas.reference         = 34.89
    wing.areas.wetted            = 57.25 
    
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[26.944,0,1.54]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 61.485 * Units.degrees  
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 0.2962
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.45
    segment.dihedral_outboard             = 0. * Units.degrees
    segment.sweeps.quarter_chord          = 31.2 * Units.degrees   
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_2'
    segment.percent_span_location         = 1.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.1183 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        

    # add to vehicle
    vehicle.append_component(wing)

    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Fuselage ############################################################## 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.number_coach_seats                 = vehicle.passengers 
    fuselage.seats_abreast                      = 6
    fuselage.seat_pitch                         = 1     * Units.meter 
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2. 
    fuselage.lengths.nose                       = 6.4   * Units.meter
    fuselage.lengths.tail                       = 8.0   * Units.meter
    fuselage.lengths.total                      = 38.02 * Units.meter  
    fuselage.lengths.fore_space                 = 6.    * Units.meter
    fuselage.lengths.aft_space                  = 5.    * Units.meter
    fuselage.width                              = 3.74  * Units.meter
    fuselage.heights.maximum                    = 3.74  * Units.meter
    fuselage.effective_diameter                 = 3.74     * Units.meter
    fuselage.areas.side_projected               = 142.1948 * Units['meters**2'] 
    fuselage.areas.wetted                       = 446.718  * Units['meters**2'] 
    fuselage.areas.front_projected              = 12.57    * Units['meters**2']  
    fuselage.differential_pressure              = 5.0e4 * Units.pascal 
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_0'    
    segment.percent_x_location                  = 0.0000
    segment.percent_z_location                  = -0.00144 
    segment.height                              = 0.0100 
    segment.width                               = 0.0100  
    fuselage.Segments.append(segment)   
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_1'    
    segment.percent_x_location                  = 0.00576 
    segment.percent_z_location                  = -0.00144 
    segment.height                              = 0.7500
    segment.width                               = 0.6500
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'   
    segment.percent_x_location                  = 0.02017 
    segment.percent_z_location                  = 0.00000 
    segment.height                              = 1.52783 
    segment.width                               = 1.20043 
    fuselage.Segments.append(segment)      
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  = 0.03170 
    segment.percent_z_location                  = 0.00000 
    segment.height                              = 1.96435 
    segment.width                               = 1.52783 
    fuselage.Segments.append(segment)   

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 0.04899 	
    segment.percent_z_location                  = 0.00431 
    segment.height                              = 2.72826 
    segment.width                               = 1.96435 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 0.07781 
    segment.percent_z_location                  = 0.00861 
    segment.height                              = 3.49217 
    segment.width                               = 2.61913 
    fuselage.Segments.append(segment)     
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 0.10375 
    segment.percent_z_location                  = 0.01005 
    segment.height                              = 3.70130 
    segment.width                               = 3.05565 
    fuselage.Segments.append(segment)             
     
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 0.16427 
    segment.percent_z_location                  = 0.01148 
    segment.height                              = 3.92870 
    segment.width                               = 3.71043 
    fuselage.Segments.append(segment)    
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'   
    segment.percent_x_location                  = 0.22478 
    segment.percent_z_location                  = 0.01148 
    segment.height                              = 3.92870 
    segment.width                               = 3.92870 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'     
    segment.percent_x_location                  = 0.69164 
    segment.percent_z_location                  = 0.01292
    segment.height                              = 3.81957
    segment.width                               = 3.81957
    fuselage.Segments.append(segment)     
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'     
    segment.percent_x_location                  = 0.71758 
    segment.percent_z_location                  = 0.01292
    segment.height                              = 3.81957
    segment.width                               = 3.81957
    fuselage.Segments.append(segment)   
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'     
    segment.percent_x_location                  = 0.78098 
    segment.percent_z_location                  = 0.01722
    segment.height                              = 3.49217
    segment.width                               = 3.71043
    fuselage.Segments.append(segment)    
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_12'     
    segment.percent_x_location                  = 0.85303
    segment.percent_z_location                  = 0.02296
    segment.height                              = 3.05565
    segment.width                               = 3.16478
    fuselage.Segments.append(segment)             
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_13'     
    segment.percent_x_location                  = 0.91931 
    segment.percent_z_location                  = 0.03157
    segment.height                              = 2.40087
    segment.width                               = 1.96435
    fuselage.Segments.append(segment)               
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_14'     
    segment.percent_x_location                  = 1.00 
    segment.percent_z_location                  = 0.04593
    segment.height                              = 1.09130
    segment.width                               = 0.21826
    fuselage.Segments.append(segment)       
    
    # add to vehicle
    vehicle.append_component(fuselage)
     

    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################################### Energy Network ##############################################################    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    #initialize the fuel network
    net                                         = RCAIDE.Framework.Networks.Fuel() 
    
    #------------------------------------------------------------------------------------------------------------------------- 
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------- 
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()  
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'starboard_propulsor'
    turbofan.active_fuel_tanks                  = ['fuel_tank']   
    turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 5.4    
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 35000.0* Units.N 
             
    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        
                   
    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 
    ram                                         = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 
          
    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle 

    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 1.9   
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.0    
    turbofan.high_pressure_compressor              = high_pressure_compressor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turbofan.low_pressure_turbine                  = low_pressure_turbine
   
    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turbofan.high_pressure_turbine                 = high_pressure_turbine 

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1500
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A1()  
    turbofan.combustor                             = combustor

    # core nozzle
    core_nozzle                                    = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turbofan.core_nozzle                           = core_nozzle
             
    # fan nozzle             
    fan_nozzle                                     = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                                 = 'fan nozzle'
    fan_nozzle.polytropic_efficiency               = 0.95
    fan_nozzle.pressure_ratio                      = 0.99 
    turbofan.fan_nozzle                            = fan_nozzle 
    
    # design turbofan
    design_turbofan(turbofan)  
    # append propulsor to distribution line  
   
 
    # Nacelle 
    nacelle                                     = RCAIDE.Library.Components.Nacelles.Body_of_Revolution_Nacelle()
    nacelle.diameter                            = 2.05
    nacelle.length                              = 2.71
    nacelle.tag                                 = 'nacelle_1'
    nacelle.inlet_diameter                      = 2.0
    nacelle.origin                              = [[13.5,4.38,-1.5]] 
    nacelle.areas.wetted                        = 1.1*np.pi*nacelle.diameter*nacelle.length 
    nacelle_airfoil                             = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    nacelle_airfoil.NACA_4_Series_code          = '2410'
    nacelle.append_airfoil(nacelle_airfoil)  
    turbofan.nacelle                            = nacelle
    
    fuel_line.propulsors.append(turbofan)  

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------      
    # copy turbofan
    turbofan_2                                  = deepcopy(turbofan)
    turbofan_2.active_fuel_tanks                = ['fuel_tank'] 
    turbofan_2.tag                              = 'port_propulsor' 
    turbofan_2.origin                           = [[13.72,-4.38,-1.1]]  # change origin 
    turbofan_2.nacelle.origin                   = [[13.5,-4.38,-1.5]]
         
    # append propulsor to distribution line 
    fuel_line.propulsors.append(turbofan_2)
  
    #------------------------------------------------------------------------------------------------------------------------- 
    #  Energy Source: Fuel Tank
    #------------------------------------------------------------------------------------------------------------------------- 
    # fuel tank
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin 
    
    # append fuel 
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Jet_A1()   
    fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel            
    
    # apend fuel tank to dataclass of fuel tanks on fuel line 
    fuel_line.fuel_tanks.append(fuel_tank) 

    # Append fuel line to Network      
    net.fuel_lines.append(fuel_line)   

    # Append energy network to aircraft 
    vehicle.append_energy_network(net)    
     
    return vehicle
  
def Concorde_vehicle_setup():
 
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Concorde'
    
    # mass properties
    vehicle.mass_properties.max_takeoff               = 185000.   # kg
    vehicle.mass_properties.operating_empty           = 78700.   # kg
    vehicle.mass_properties.takeoff                   = 183000.   # kg, adjusted due to significant fuel burn on runway
    vehicle.mass_properties.cargo                     = 1000.  * Units.kilogram   
    vehicle.mass_properties.max_zero_fuel             = 92000.
        
    # envelope properties
    vehicle.envelope.ultimate_load = 3.75
    vehicle.envelope.limit_load    = 2.5

    # basic parameters
    vehicle.reference_area               = 358.25      
    vehicle.passengers                   = 100
    vehicle.systems.control              = "fully powered" 
    vehicle.systems.accessories          = "sst"
    vehicle.maximum_cross_sectional_area = 13.9
    vehicle.total_length                 = 61.66
    vehicle.design_mach_number           = 2.02
    vehicle.design_range                 = 4505 * Units.miles
    vehicle.design_cruise_alt            = 60000.0 * Units.ft
    
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
    rel_path                       = os.path.dirname(ospath) + separator + '..' + separator    
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
    fuel_line.propulsors.append(outer_right_turbojet) 

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
    fuel_line.propulsors.append(inner_right_turbojet) 

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
    fuel_line.propulsors.append(inner_left_turbojet) 

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
    fuel_line.propulsors.append(outer_left_turbojet) 
 
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
    
     # Append fuel line to network      
    net.fuel_lines.append(fuel_line)    
  
    #------------------------------------------------------------------------------------------------------------------------------------          
    # Append energy network to aircraft 
    vehicle.append_energy_network(net)     


    return vehicle 
 
if __name__ == '__main__': 
    main()    
    plt.show()