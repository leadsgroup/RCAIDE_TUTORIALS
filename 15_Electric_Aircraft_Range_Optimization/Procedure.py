# Procedure.py
# 
# Created:  Mar 2023, M Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------      
from RCAIDE.Analyses.Process import Process  
import  Missions

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------       

def setup():
    procedure = Process() 
    procedure.modify_vehicle_and_mission  = modify_vehicle_and_mission
    procedure.finalize                    = finalize 
    procedure.missions                    = Process()
    procedure.missions.base               = evaluate_mission 
    procedure.post_process                = post_process
        
    return procedure 

# ----------------------------------------------------------------------        
#   Evaluate Mission
# ----------------------------------------------------------------------    
    
def evaluate_mission(nexus): 
    # Evaluate the missions and save to results    
    mission         = nexus.missions.base
    results         = nexus.results
    results.base   = mission.evaluate()                                        
    return nexus
 
# ----------------------------------------------------------------------        
#  Modify Vehicle
# ----------------------------------------------------------------------    

def modify_vehicle_and_mission(nexus): 
    # this function actually does nothing but here is where you would modify the vehicle/mission
    # segments that rely on the mission  
    # Pull out the vehicle
    vehicle_opt = nexus.vehicle_configurations.base
 
    # ----------------------------------------------------------------------
    # Update Mission 
    # ---------------------------------------------------------------------- 
    # Re-set nattery charge each optimization
    nexus.missions.base.segments.takeoff.battery_energy =  vehicle_opt.networks.battery_electric_rotor.battery.pack.max_energy 

    # diff the new data
    vehicle_opt.store_diff()
    
    return nexus 
 

# ----------------------------------------------------------------------
#   Finalizing Function
# ----------------------------------------------------------------------    

def finalize(nexus):
    
    nexus.analyses.finalize()   
    
    return nexus         

# ----------------------------------------------------------------------
#   Post Process Results to give back to the optimizer
# ----------------------------------------------------------------------   

def post_process(nexus):
    
    # Unpack data  
    res       = nexus.results.base
    final_SOC = res.segments[-1].conditions.propulsion.battery.cell.state_of_charge[-1,0]
    
    summary         = nexus.summary     
    summary.Nothing = 0.0 
    summary.SOC_EOF = final_SOC           
  
    nexus.total_number_of_iterations +=1
    return nexus 