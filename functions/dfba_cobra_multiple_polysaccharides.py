import cobra
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from collections import OrderedDict

import copy
from sklearn.metrics import mean_squared_error 

name2id = {"Glucose":"glc__D_e",
                               "Xylose":"xyl__D_e",
                               "Galactose":"gal_e",
                               "Arabinose":"arab__L_e",
                              "Glucuronic acid":"glcur_c"}


def add_dynamic_bounds(model, conc_dict,glc_eq_poly_dict,vmax_inner_glc,Km_inner_glc,vmax_outer,Km_outer):
    """Use external concentrations to bound the uptake flux of glucose."""    
    medium = model.medium
    
    cellulase_dict = OrderedDict()
    for polysac_id,glc_eq_dict in  glc_eq_poly_dict.items():
        
        for met_id,glc_eq in  glc_eq_dict.items():
            
            max_import = -vmax_inner_glc*conc_dict[met_id]/ ((Km_inner_glc + conc_dict[met_id]))#*glc_eq) 
                
            medium[met_id]=-max_import

        cellulase = -vmax_outer*conc_dict[polysac_id]/(Km_outer + conc_dict[polysac_id])
    
        cellulase_real=cellulase

        cellulase_dict[polysac_id] = cellulase_real
        

    model.medium=medium
    return cellulase_dict


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


def dynamic_system_general(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    """Calculate the time derivative of external species."""    

    vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer = combination

    rxns_map = copy.copy(rxns)
    
    for key in glc_eq_poly_dict.keys():
        rxns_map.remove(key)
        
    conc_dict = dict(zip(rxns,y))
    # Calculate the specific exchanges fluxes at the given external concentrations.

    with model:
        # Calculate the specific exchanges fluxes at the given external concentrations.
        cellulase_dict = add_dynamic_bounds(model, conc_dict,glc_eq_poly_dict,vmax_inner_glc, Km_inner_glc,vmax_outer, Km_outer)
        feasibility = cobra.util.fix_objective_as_constraint(model)
        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns_map, objective_dir)

    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.

    fluxes =lex_constraints.values
    i = 1 

    for polysac_id,glc_eq_dict in  glc_eq_poly_dict.items():  
        for key, glc_eq in glc_eq_dict.items():
            uptake = fluxes[i]
            fluxes[i] = uptake - cellulase_dict[polysac_id]/(glc_eq*len(glc_eq_dict)) # 1 mol of cellulase produces one mol of glucose equivalents
            i +=1
        fluxes =np.append(fluxes, cellulase_dict[polysac_id])

    fluxes *= conc_dict["Growth"]
    return fluxes


def multiple_polysaccharide_inner_problem(combination,model,media,rxns,y0,objective_dir,glc_eq_poly_dict,t_end=False):
    vmax_inner_glc, Km_inner_glc,vmax_outer, Km_outer= combination
    
    y = []
    
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
        events = [infeasible_event,infeasible_event2,infeasible_event3,infeasible_event4,infeasible_event5],
        args = (model,rxns,objective_dir,glc_eq_poly_dict,combination),
    )

    return sol




def multiple_polysaccharide_simulation(combination,model,media,rxns,y0,objective_dir,glc_eq_poly_dict,t_end=False):

    
    C_dict = OrderedDict(zip(rxns,y0))
    C_results_tot = {rxn:{"C":[value],"t":[0]} for rxn,value in C_dict.items()}
    
    glc_eq_poly_dict_copy = glc_eq_poly_dict.copy()
    C_dict_copy = C_dict.copy()
    objective_dir_copy = objective_dir.copy()
    
    


    while sum([C_results_tot[polysac_id]["C"][-1] for polysac_id in glc_eq_poly_dict_copy.keys()]) >1e-1:
        
        sol = multiple_polysaccharide_inner_problem(combination,
                                                        model,
                                                        media,
                                                        rxns,
                                                        y0,
                                                        objective_dir_copy,
                                                        glc_eq_poly_dict_copy,
                                                        t_end=t_end)
        
            
        ## Update concentrations
        C_results = dict(zip(rxns,sol.y))

        for rxn,conc in C_results.items():
            C_results_tot[rxn]["C"] = C_results_tot[rxn]["C"][:-1] + list(conc)
            C_results_tot[rxn]["t"] = C_results_tot[rxn]["t"][:-1] + list(C_results_tot[rxn]["t"][-1] +sol.t)


        last_events = [(t_event[-1] if len(t_event) else 0) for t_event in sol.t_events]
        event_nr = np.argmax(last_events)

        ## Prepare for next iteration 
        if event_nr == 1:
            polysac_id = "EX_arabinoxylan_e"

        elif event_nr== 2:
            polysac_id = "EX_xyloglucan_e"   
        elif event_nr== 3:
            polysac_id = "EX_xylan_e"   

        elif event_nr== 4:
            polysac_id = "EX_cellulose_e"   
        else:
            break

        # remove items that have reached the final concentration
        rm_rxns_dict = glc_eq_poly_dict_copy.pop(polysac_id)
        rm_rxns = [polysac_id]+list(rm_rxns_dict.keys())

        for rxn in rm_rxns: C_dict_copy.pop(rxn); objective_dir_copy.pop()

        objective_dir_copy.append('max')
        rxns = list(C_dict_copy.keys())
        y0 = np.array([C_results[key][-1] for key in C_dict_copy.keys()])

    return C_results_tot



def infeasible_event(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer = combination
    rxns_map = copy.copy(rxns)
    
    for key in glc_eq_poly_dict.keys():
        rxns_map.remove(key)
    conc_dict = dict(zip(rxns,y))
     
    with model:

        cellulase_dict = add_dynamic_bounds(model, conc_dict,glc_eq_poly_dict,vmax_inner_glc, Km_inner_glc, vmax_outer, Km_outer)
        feasibility = cobra.util.fix_objective_as_constraint(model)
    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True

def infeasible_event2(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    
    conc_dict = dict(zip(rxns,y))
    
    if "EX_arabinoxylan_e" in rxns:
        arabinoxylan = conc_dict["EX_arabinoxylan_e"]+conc_dict["EX_AX_e"]+conc_dict["EX_AXX_e"]+conc_dict["EX_XAXX_e"]+conc_dict["EX_A23XX_e"] + +conc_dict["EX_XA23XX_e"]
        
        diff = arabinoxylan - infeasible_event2.epsilon
        if diff<0:
            print(f"Lack of arabinoxylan: {diff}")
        return diff
    else:
        return 10

infeasible_event2.epsilon=1E-2
infeasible_event2.direction = 0
infeasible_event2.terminal = True

def infeasible_event3(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    
    conc_dict = dict(zip(rxns,y))
    
    if "EX_xyloglucan_e" in rxns: 
        xyloglucan = conc_dict["EX_xyloglucan_e"]+conc_dict["EX_QLLG_e"]+conc_dict["EX_QLQG_e"]+conc_dict["EX_QQLG_e"]+conc_dict["EX_QQQG_e"]+conc_dict["EX_GQQG_e"]
        diff = xyloglucan - infeasible_event3.epsilon
        if diff<0:
            print(f"Lack of xyloglucan: {diff}")
        return diff
    else:
        return 10

infeasible_event3.epsilon=1E-2
infeasible_event3.direction = 0
infeasible_event3.terminal = True

def infeasible_event4(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    
    conc_dict = dict(zip(rxns,y))
    
    if "EX_xylan_e" in rxns: 
        xylan = conc_dict["EX_xylan_e"]+conc_dict["EX_xylb_e"]+conc_dict["EX_xyl3_e"]+conc_dict["EX_xylan4_e"]+conc_dict["EX_xylan8_e"]
        diff = xylan - infeasible_event4.epsilon
        if diff<0:
            print(f"Lack of xylan: {diff}")
        return diff
    else:
        return 10

infeasible_event4.epsilon=1E-2
infeasible_event4.direction = 0
infeasible_event4.terminal = True

def infeasible_event5(t, y,model,rxns,objective_dir,glc_eq_poly_dict,combination):
    
    conc_dict = dict(zip(rxns,y))
    
    if "EX_cellulose_e" in rxns: 
        cellulose = conc_dict["EX_cellulose_e"]+conc_dict["EX_cell5_e"]+conc_dict["EX_cell4_e"]+conc_dict["EX_cell3_e"]+conc_dict["EX_cellb_e"]
        diff = cellulose - infeasible_event5.epsilon
        if diff<0:
            print(f"Lack of cellulose: {diff}")
        return diff
    else:
        return 10

infeasible_event5.epsilon=1E-2
infeasible_event5.direction = 0
infeasible_event5.terminal = True




