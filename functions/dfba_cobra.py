import cobra
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

import copy
from sklearn.metrics import mean_squared_error 



def add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc,Km_inner_glc,vmax_outer,Km_outer):
    """Use external concentrations to bound the uptake flux of glucose."""    
    medium = model.medium
    
    for met_id,glc_eq in  glc_eq_dict.items():
        max_import = -vmax_inner_glc*conc_dict[met_id]/ ((Km_inner_glc/glc_eq + conc_dict[met_id])*glc_eq) ## 
        medium[met_id]=-max_import
        

    cellulase = -vmax_outer*conc_dict["EX_polysac_e"]/(Km_outer + conc_dict["EX_polysac_e"])
    
    #print(f"cellulase: {-cellulase}, glc_uptake: {medium['EX_glc__D_e']}, cellb_uptake: {medium['EX_cellb_e']}")
    model.medium=medium
    return cellulase


def read_model(media,lp_feasibility=True):
    model = cobra.io.read_sbml_model('../models/RcH10_final.xml')
    exchanges_ids = [rxn.id for rxn in model.exchanges]
    media_exchanges = ["EX_"+met+"_e" for met in media["DM_cellobiose"]] 
    medium = model.medium

    for exchange in exchanges_ids:
        if exchange in ["EX_xyl__D_e","EX_glc__D_e","EX_cellb_e","EX_gal_e","EX_gal__L_e"]:
            medium[exchange]=0

        elif exchange in media_exchanges:
            medium[exchange]=100 
        else: 
            medium[exchange]=0

    model.medium = medium
    if lp_feasibility:
        cobra.util.add_lp_feasibility(model)
    return model


def dynamic_system_general(t, y,model,rxns,objective_dir,glc_eq_dict,combination):
    """Calculate the time derivative of external species."""    

    vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer = combination

    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_polysac_e")
    conc_dict = dict(zip(rxns,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        # Calculate the specific exchanges fluxes at the given external concentrations.
        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc, Km_inner_glc,vmax_outer, Km_outer)
        feasibility = cobra.util.fix_objective_as_constraint(model)
        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns_map, objective_dir)
        
        
    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    
    
    fluxes =lex_constraints.values
    #fluxes =[0,0,0]
    i = 1 
    
    for key, glc_eq in glc_eq_dict.items():
        uptake = fluxes[i]
        fluxes[i] = uptake - cellulase/(glc_eq*len(glc_eq_dict)) # 1 mol of cellulase produces one mol of glucose equivalents
        i +=1
        
    #    if key=="EX_glc__D_e": print(f"flux glc: {fluxes[-1]}")
    #    if key=="EX_cellb_e": print(f"flux cellb: {fluxes[-1]}")
        
    fluxes =np.append(fluxes, cellulase)
    fluxes *= conc_dict["Growth"]

    return fluxes



def optimize_parameters_inner_problem_general(combination,model,media,rxns,y0,objective_dir,glc_eq_dict,alternative_solution=False,t_end=False):
    print(combination)

    vmax_inner_glc, Km_inner_glc,vmax_outer, Km_outer= combination
    
    #try:
    if t_end:
        ts = np.linspace(0, t_end, 1000)   
    else:
        ts = np.linspace(0, 50, 1000)   

    sol = solve_ivp(
        fun=dynamic_system_general,
        t_span=(ts.min(), ts.max()),
        y0=y0,
        t_eval=ts,
        method='LSODA',
        events = [infeasible_event],
        args = (model,rxns,objective_dir,glc_eq_dict,combination)
    )
    #except:
    #    print(f"had issues with this combination: {combination}")

    return sol
    

def infeasible_event(t, y,model,rxns,objective_dir,glc_eq_dict,combination):
    vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer = combination
    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_polysac_e")
    conc_dict = dict(zip(rxns,y))
     
    with model:

        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer)
        feasibility = cobra.util.fix_objective_as_constraint(model)
    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True
