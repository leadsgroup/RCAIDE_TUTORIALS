# Analyses.py
# 
# Created:  Mar 2016, M. Vegh
# Modified: Aug 2017, E. Botero

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import Units

import numpy as np

# ----------------------------------------------------------------------        
#   Setup Analyses
# ----------------------------------------------------------------------  

def setup(configs):
    
    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base(config)
        if tag == 'cruise_spoilers':
            # this is done since drag is not sufficient for the desired profile
            analysis.aerodynamics.settings.spoiler_drag_increment = 0.005       
        analyses[tag] = analysis


    return analyses

# ----------------------------------------------------------------------        
#   Define Base Analysis
# ----------------------------------------------------------------------  

def base(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = RCAIDE.Framework.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.number_of_spanwise_vortices  = 5 
    aerodynamics.settings.number_of_chordwise_vortices = 1
    aerodynamics.process.compute.lift.inviscid_wings.training.angle_of_attack = np.array([[-5., 0.0, 5.0, 10.0, 75.]]).T * Units.deg 
    aerodynamics.process.compute.lift.inviscid_wings.training.Mach            = np.array([[0.0, 0.2, 0.5, 0.70, 0.80, 0.9, 1.3, 1.35, 1.5, 2.0]]).T         
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()
    stability.vehicle = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Framework.Analyses.Energy.Energy()
    energy.networks = vehicle.networks
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

    return analyses 