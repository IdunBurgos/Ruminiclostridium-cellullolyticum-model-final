import cobra
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

import copy
from sklearn.metrics import mean_squared_error 



cellulose_exp = pd.read_csv("../input/Desvaux2001_batch_data/cellulose_g.csv")
molar_mass = 180.156  # Based on glucose equivalent (based on 0.35 cellobiose and 0.3 glucose)
cellulose_exp["y mmol"]= cellulose_exp[" y"].map(lambda x: (x/molar_mass)*1000)
cellulose_exp = cellulose_exp[cellulose_exp.x<70].copy()


biomass_exp = pd.read_csv("../input/Desvaux2001_batch_data/biomass_mg.csv")
biomass_exp["y g"] = biomass_exp[" y"].map(lambda x: x/1000)
biomass_exp = biomass_exp[biomass_exp.x<70].copy()

cellb_exp = pd.read_csv("../input/Desvaux2001_batch_data/cellobiose_mmol.csv")
cellb_exp = cellb_exp[cellb_exp.x<70].copy()

glc_exp = pd.read_csv("../input/Desvaux2001_batch_data/glucose_mmol.csv")
glc_exp = glc_exp[glc_exp.x<70].copy()

acetate_exp = pd.read_csv("../input/Desvaux2001_batch_data/acetate_mmol.csv")
acetate_exp = acetate_exp[acetate_exp.x<70].copy()

ethanol_exp = pd.read_csv("../input/Desvaux2001_batch_data/ethanol_mmol.csv")
ethanol_exp = ethanol_exp[ethanol_exp.x<70].copy()

lactate_exp = pd.read_csv("../input/Desvaux2001_batch_data/lactate_mmol.csv")
lactate_exp = lactate_exp[lactate_exp.x<70].copy()

def add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc,Km_inner_glc,vmax_inner_cellb,Km_inner_cellb,vmax_outer,Km_outer,Ki):
    """Use external concentrations to bound the uptake flux of glucose."""    
    medium = model.medium
    

    max_import = -vmax_inner_glc*conc_dict["EX_glc__D_e"]/ (Km_inner_glc + conc_dict["EX_glc__D_e"]) ## 
    medium["EX_glc__D_e"]=-max_import
    
    max_import = -vmax_inner_cellb*conc_dict["EX_cellb_e"]/ (Km_inner_cellb + conc_dict["EX_cellb_e"]) ## 
    medium["EX_cellb_e"]=-max_import
        
    cellulase = -vmax_outer*(conc_dict["EX_cellulose_e"]/((1 + (conc_dict["EX_cellb_e"]/Ki))*Km_outer + conc_dict["EX_cellulose_e"]))
    
    model.medium=medium
    return cellulase


def read_model(media,lp_feasibility=True):
    model = cobra.io.read_sbml_model('../models/RcH10_final_flux_ratio.xml')
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


def dynamic_system_cellulose(t, y,model,rxns,objective_dir,glc_eq_dict,combination):
    """Calculate the time derivative of external species."""    

    vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki = combination

    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_cellulose_e")
    conc_dict = dict(zip(rxns,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        # Calculate the specific exchanges fluxes at the given external concentrations.
        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki)
        feasibility = cobra.util.fix_objective_as_constraint(model)
        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns_map, objective_dir)
        
    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    fluxes =lex_constraints.values
    i = 1 
    for key, glc_eq in glc_eq_dict.items():
        uptake = fluxes[i]
        fluxes[i] = uptake - cellulase/(glc_eq*len(glc_eq_dict)) # 1 mol of cellulase produces one mol of glucose equivalents
        i +=1
        
    fluxes =np.append(fluxes, cellulase)
    fluxes *= conc_dict["Growth"]

    return fluxes



def optimize_parameters_inner_problem_cellulose(combination,model,media,rxns,y0,objective_dir,glc_eq_dict,alternative_solution=False,t_end=False):
    print(combination)

    vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki= combination
    
    try:
        if t_end:
            ts = np.linspace(biomass_exp.iloc[0,0], t_end, 1000)   
        else:
            ts = np.linspace(biomass_exp.iloc[0,0], biomass_exp.iloc[biomass_exp["x"].size-1,0], 1000)   

        sol = solve_ivp(
            fun=dynamic_system_cellulose,
            t_span=(ts.min(), ts.max()),
            y0=y0,
            t_eval=ts,
            method='LSODA',
            events = [infeasible_event],
            args = (model,rxns,objective_dir,glc_eq_dict,combination)
        )
    except:
        print(f"had issues with this combination: {combination}")

    
    C_dict_results = dict(zip(rxns,sol.y))
    
    # Growth
    
    interp_func = interp1d(sol.t,C_dict_results["Growth"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(biomass_exp.x.values)
    
    score_spec = mean_squared_error(biomass_exp["y g"].values,y_interp)/(biomass_exp["y g"].mean()) 
    penalty_growth =1000*(score_spec)
    
    
    #  Cellulose
    interp_func = interp1d(sol.t, C_dict_results["EX_cellulose_e"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(cellulose_exp.x)
    
    score_spec = mean_squared_error(cellulose_exp["y mmol"],y_interp)/(cellulose_exp["y mmol"].mean()) 
    penalty_cellulose =100*(score_spec)
    
    #  Acetate    
    interp_func = interp1d(sol.t, C_dict_results["EX_ac_e"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(acetate_exp.x)
    
    score_spec = mean_squared_error(acetate_exp[" y"],y_interp)/(acetate_exp[" y"].mean()) 
    penalty_acetate =10*(score_spec)

    #  Ethanol    
    interp_func = interp1d(sol.t, C_dict_results["EX_etoh_e"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(ethanol_exp.x)
    
    score_spec = mean_squared_error(ethanol_exp[" y"],y_interp)/(ethanol_exp[" y"].mean()) 
    penalty_etoh =10*(score_spec)
    
    
    #  Lactate    
    interp_func = interp1d(sol.t, C_dict_results["EX_lac__L_e"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(lactate_exp.x)
    
    score_spec = mean_squared_error(lactate_exp[" y"],y_interp)/(lactate_exp[" y"].mean()) 
    penalty_lac =10*(score_spec)
        
    # Penalty
        
    penalty = penalty_growth + penalty_cellulose + penalty_acetate + penalty_etoh + penalty_lac
    
    print(f"penalty: {penalty}")
    
    if alternative_solution:
        return sol,penalty,{"penalty_growth":penalty_growth,"penalty_cellulose":penalty_cellulose,"penalty_acetate":penalty_acetate,"penalty_ethanol":penalty_etoh,"penalty_lac":penalty_lac}
    else:
        return penalty
    
    

def infeasible_event(t, y,model,rxns,objective_dir,glc_eq_dict,combination):

    
    #vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki,ratio_ac_etoh,ratio_ac_lac = combination
    vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki = combination

    
    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_cellulose_e")
    conc_dict = dict(zip(rxns,y))
     
    
    with model:

        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb, vmax_outer, Km_outer,Ki)
        feasibility = cobra.util.fix_objective_as_constraint(model)
    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True
