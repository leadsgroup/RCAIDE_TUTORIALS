'''

The script below documents how to set up and plot the results of an isolaed/static propeller analysis  

''' 
#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Framework.Core                              import Units
from RCAIDE.Library.Plots                               import *    
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor import design_propeller 
from RCAIDE.Framework.Mission.Common                    import Results  
from RCAIDE.Framework.Mission.Segments.Segment          import Segment 
from RCAIDE.Framework.Mission.Common                    import Conditions
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor.compute_rotor_performance import compute_rotor_performance 

import os
import numpy as np 
import matplotlib.pyplot as plt   

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main(): 
    # design aircract 
    prop = design_test_propeller()
    
    # operarting states 
    fligth_path_angle =  0.0
    true_course       =  0.0
    ctrl_pts          =  1
    AoA               =  0.0 
    
    # Find the operating conditions
    atmosphere                                             = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions                                  = atmosphere.compute_values(prop.cruise.design_altitude)  
 
    segment                                                = Segment()  
    conditions                                             = Results()  
    conditions.aerodynamics.angle_of_attack                = np.atleast_2d(AoA).T
    conditions.freestream.density                          = np.ones((ctrl_pts,1)) * atmosphere_conditions.density[0][0]
    conditions.freestream.dynamic_viscosity                = np.ones((ctrl_pts,1)) * atmosphere_conditions.dynamic_viscosity[0][0]   
    conditions.freestream.speed_of_sound                   = np.ones((ctrl_pts,1)) * atmosphere_conditions.speed_of_sound[0][0]
    conditions.freestream.temperature                      = np.ones((ctrl_pts,1)) * atmosphere_conditions.temperature[0][0]  
    conditions.frames.inertial.velocity_vector             = np.array([[prop.cruise.design_freestream_velocity,0,0]]) 
    conditions.frames.planet.true_course                   = np.zeros((ctrl_pts,3,3)) 
    conditions.frames.planet.true_course[:,0,0]            = np.cos(true_course),
    conditions.frames.planet.true_course[:,0,1]            = - np.sin(true_course)
    conditions.frames.planet.true_course[:,1,0]            = np.sin(true_course)
    conditions.frames.planet.true_course[:,1,1]            = np.cos(true_course) 
    conditions.frames.planet.true_course[:,2,2]            = 1 
    conditions.frames.wind.transform_to_inertial           = np.zeros((ctrl_pts,3,3))   
    conditions.frames.wind.transform_to_inertial[:,0,0]    = np.cos(fligth_path_angle) 
    conditions.frames.wind.transform_to_inertial[:,0,2]    = np.sin(fligth_path_angle) 
    conditions.frames.wind.transform_to_inertial[:,1,1]    = 1 
    conditions.frames.wind.transform_to_inertial[:,2,0]    = -np.sin(fligth_path_angle) 
    conditions.frames.wind.transform_to_inertial[:,2,2]    = np.cos(fligth_path_angle)  
    conditions.frames.body.transform_to_inertial           = np.zeros((ctrl_pts,3,3))
    conditions.frames.body.transform_to_inertial[:,0,0]    = np.cos(AoA)
    conditions.frames.body.transform_to_inertial[:,0,2]    = np.sin(AoA)
    conditions.frames.body.transform_to_inertial[:,1,1]    = 1
    conditions.frames.body.transform_to_inertial[:,2,0]    = -np.sin(AoA)
    conditions.frames.body.transform_to_inertial[:,2,2]    = np.cos(AoA)  
    segment.state.conditions                               = conditions 

    bus                                      = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus() 
    electric_rotor                           = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    electric_rotor.rotor                     = prop
    bus.propulsors.append(electric_rotor)    
    segment.state.conditions.energy[bus.tag] = Conditions()
    segment.state.conditions.noise[bus.tag]  = Conditions()
    electric_rotor.append_operating_conditions(segment,bus) 
    for tag, item in  electric_rotor.items(): 
        if issubclass(type(item), RCAIDE.Library.Components.Component):
            item.append_operating_conditions(segment,bus,electric_rotor)
            
    # Run BEMT
    segment.state.conditions.expand_rows(ctrl_pts)
    rotor_conditions             =  segment.state.conditions.energy[bus.tag][electric_rotor.tag][prop.tag]     
    rotor_conditions.omega[:,0]  = prop.cruise.design_angular_velocity
    compute_rotor_performance(electric_rotor,segment.state,bus)
    
    plot_results(conditions.energy[bus.tag][electric_rotor.tag][prop.tag], prop,'blue','-','s')
 
    return

def design_test_propeller(): 
    
    prop                                     = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    prop.number_of_blades                    = 3
    prop.number_of_engines                   = 1
    prop.tip_radius                          = 1.0668
    prop.hub_radius                          = 0.21336
    prop.cruise.design_freestream_velocity   = 49.1744
    prop.cruise.design_tip_mach              = 0.65
    prop.cruise.design_angular_velocity      = 207.16160479940007
    prop.cruise.design_Cl                    = 0.7
    prop.cruise.design_altitude              = 1. * Units.km 
    prop.cruise.design_thrust                = 3054.4809132125697
    
    # define first airfoil    
    ospath                                     = os.path.abspath(__file__)
    separator                                  = os.path.sep
    rel_path                                   = os.path.dirname(ospath) + separator + '..'+ separator    
    airfoil_1                                  = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil_1.tag                              = 'NACA_4412' 
    airfoil_1.coordinate_file                  = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil_1.polar_files                      =[rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt'] 
    prop.append_airfoil(airfoil_1)           # append first airfoil 
    
    # define  second airfoil 
    airfoil_2                                = RCAIDE.Library.Components.Airfoils.Airfoil()  
    airfoil_2.tag                            = 'Clark_Y' 
    airfoil_2.coordinate_file                =   rel_path + 'Airfoils' + separator + 'Clark_y.txt' 
    airfoil_2.polar_files                    = [ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'Clark_y_polar_Re_50000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'Clark_y_polar_Re_100000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'Clark_y_polar_Re_200000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'Clark_y_polar_Re_500000.txt',
                                                 rel_path + 'Airfoils' + separator + 'Polars' + separator + 'Clark_y_polar_Re_1000000.txt'] 
    prop.append_airfoil(airfoil_2)          # append second airfoil 
    
    # define polar stations on rotor 
    prop.airfoil_polar_stations           = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
    design_propeller(prop)
    
    # plot propeller 
    plot_3d_rotor(prop) 
    
    return prop 
     
def plot_results(results,prop,c,ls,m):

    tag                = prop.tag
    va_ind             = results.blade_axial_induced_velocity[0]  
    vt_ind             = results.blade_tangential_induced_velocity[0]  
    r                  = prop.radius_distribution
    T_distribution     = results.blade_thrust_distribution[0] 
    vt                 = results.blade_tangential_velocity[0]  
    va                 = results.blade_axial_velocity[0] 
    Q_distribution     = results.blade_torque_distribution[0] 

    # ----------------------------------------------------------------------------
    # 2D - Plots  Plots    
    # ---------------------------------------------------------------------------- 
    # perpendicular velocity, up Plot 
    fig = plt.figure('va_ind')         
    plt.plot(r  , va_ind ,color = c  , marker = m, linestyle = ls , label =  tag)          
    plt.xlabel('Radial Location')
    plt.ylabel('Induced Axial Velocity') 
    plt.legend(loc='lower right') 

    fig = plt.figure('vt_ind')          
    plt.plot(r  , vt_ind ,color = c ,marker = m, linestyle = ls , label =  tag )       
    plt.xlabel('Radial Location')
    plt.ylabel('Induced Tangential Velocity') 
    plt.legend(loc='lower right')  

    fig = plt.figure('T')     
    plt.plot(r , T_distribution ,color = c ,marker = m, linestyle = ls, label =  tag  )    
    plt.xlabel('Radial Location')
    plt.ylabel('Trust, N')
    plt.legend(loc='lower right')

    fig = plt.figure('Q')
    plt.plot(r , Q_distribution ,color = c ,marker = m, linestyle = ls, label =  tag)            
    plt.xlabel('Radial Location')
    plt.ylabel('Torque, N-m')
    plt.legend(loc='lower right')

    fig = plt.figure('Va')     
    plt.plot(r , va ,color = c  ,marker =m, linestyle = ls, label =  tag + 'axial vel')          
    plt.xlabel('Radial Location')
    plt.ylabel('Axial Velocity') 
    plt.legend(loc='lower right') 

    fig = plt.figure('Vt')       
    plt.plot(r , vt ,color = c ,marker = m, linestyle = ls, label =  tag )         
    plt.xlabel('Radial Location')
    plt.ylabel('Tangential Velocity') 
    plt.legend(loc='lower right')  

    return 

# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()
