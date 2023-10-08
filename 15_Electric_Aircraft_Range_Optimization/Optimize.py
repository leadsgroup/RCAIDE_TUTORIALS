# Optimize.py
# Created:  Mar 2023, M. Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
assert RCAIDE.__version__=='1.0.0', 'These tutorials only work with the RCAIDE 1.0.0 release'
from RCAIDE.Core import Units, Data
import numpy as np
import matplotlib.pyplot as plt
import Vehicles
import Analyses
import Missions
import Procedure
import Plot_Mission 
from RCAIDE.Optimization import Nexus 
import RCAIDE.Optimization.Package_Setups.scipy_setup as scipy_setup 

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():
    
    problem = setup()  
    output = scipy_setup.SciPy_Solve(problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-3)
    print (output)     
    
    Plot_Mission.plot_mission(problem)
    
    return 
 

def setup(): 

    nexus   = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,  (lb , ub)    , scaling , units ]  
    problem.inputs = np.array([
                  [ 'cruise_range' , 130   , 10 , 200.  ,   100. , 1*Units.nmi  ],  
    ],dtype=object)    

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([
        [ 'SOC_end_of_flight', 1. , 1*Units.less ] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
        [ 'SOC_end_of_flight' , '>' , 0.2 ,  1.  , 1*Units.less], 
         
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ]

    problem.aliases = [
        [ 'cruise_range'              , 'missions.base.segments.cruise.distance'],  
        [ 'SOC_end_of_flight'         , 'summary.SOC_EOF'],
        [ 'Nothing'                   , 'summary.Nothing'],               
    ]    
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.setup()
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.setup(nexus.analyses)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure =  Procedure.setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()    
    nexus.total_number_of_iterations = 0
    return nexus    

if __name__ == '__main__':
    main()
    plt.show()
