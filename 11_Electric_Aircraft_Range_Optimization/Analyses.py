# Analyses.py
# 
# Created:  Mar 2023, M. Clarke

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE  

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle):

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
    weights = RCAIDE.Framework.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle  
    aerodynamics.settings.model_fuselage = True 
    aerodynamics.settings.number_of_spanwise_vortices           = 25
    aerodynamics.settings.number_of_chordwise_vortices          = 5    
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