'''

The script below documents how to set up and plot the results of an isolaed/static propeller analysis  

''' 
#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Framework.Core import Units, Data 
from RCAIDE.Library.Plots  import *    
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor import design_propeller 

import os
import numpy as np 
import matplotlib.pyplot as plt   

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main(): 
 
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
    prop                                  = design_propeller(prop)

    # plot propeller 
    plot_3d_rotor(prop) 

    # Find the operating conditions
    atmosphere                                          = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions                               = atmosphere.compute_values(prop.cruise.design_altitude)  
    conditions                                          = RCAIDE.Framework.Mission.Common.Results()
    conditions._size                                    = 1
    conditions.freestream                               = Data()
    conditions.propulsion                               = Data()
    conditions.frames                                   = Data()
    conditions.frames.body                              = Data()
    conditions.frames.inertial                          = Data()
    conditions.freestream.update(atmosphere_conditions)
    conditions.freestream.dynamic_viscosity             = atmosphere_conditions.dynamic_viscosity
    conditions.frames.inertial.velocity_vector          = np.array([[prop.cruise.design_freestream_velocity,0,0]])
    conditions.propulsion.throttle                      = np.array([[1.0]])
    conditions.frames.body.transform_to_inertial        = np.array([np.eye(3)])

    conditions.frames.inertial.velocity_vector   = np.array([[prop.cruise.design_freestream_velocity,0,0]]) 

    # Create and attach this propeller 
    prop.inputs.omega  = np.array(prop.cruise.design_angular_velocity,ndmin=2) 
 
    prop.inputs.pitch_command                = 0.0*Units.degree
    F_a, Q_a, P_a, Cplast_a ,output_a , etap_a = prop.spin(conditions)  
    plot_results(output_a, prop,'blue','-','s')
 
    return
 
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
