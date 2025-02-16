import cobra
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

import copy
from sklearn.metrics import mean_squared_error 




def add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner,Km_inner,vmax_outer,Km_outer,Ki):
    """Use external concentrations to bound the uptake flux of glucose."""    
    medium = model.medium
    
    for key in glc_eq_dict.keys():
        
        max_import = -vmax_inner * conc_dict[key] / (Km_inner + conc_dict[key])
        medium[key]=-max_import
        
    if Ki==0: ### 'Simplified' 3
        cellulase = -vmax_outer*conc_dict["EX_xylogluc_e"]/(Km_outer + conc_dict["EX_xylogluc_e"])
    else:
        cellulase = -vmax_outer*conc_dict["EX_xylogluc_e"]/((1 + (conc_dict["EX_QQQG_e"]/Ki))*Km_outer + conc_dict["EX_xylogluc_e"])

    model.medium=medium
    return cellulase


def read_model(media):
    model = cobra.io.read_sbml_model('../models/RcH10_final.xml')
    exchanges_ids = [rxn.id for rxn in model.exchanges]

    media_exchanges = ["EX_"+met+"_e" for met in media["DM_xyloglucan"]]
    medium = model.medium

    for exchange in exchanges_ids:
        if exchange in ["EX_xyl__D_e","EX_glc__D_e","EX_cellb_e","EX_gal_e","EX_gal__L_e"]:
            medium[exchange]=0

        elif exchange in media_exchanges:
            medium[exchange]=100 

        else: 
            medium[exchange]=0

    model.medium = medium
    cobra.util.add_lp_feasibility(model)
    return model
    
    
def dynamic_system_oligo(t, y,model,rxns,objective_dir,glc_eq_dict,combination):
    """Calculate the time derivative of external species."""    
    vmax_inner, Km_inner = combination
    vmax_outer, Km_outer,Ki = [0.1,0.1,0.1]
    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_xylogluc_e")
    conc_dict = dict(zip(rxns,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner,Km_inner,vmax_outer,Km_outer,Ki)
        feasibility = cobra.util.fix_objective_as_constraint(model)
        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns_map, objective_dir)
        
    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    fluxes =lex_constraints.values
    i = 1 
    for key, glc_eq in glc_eq_dict.items():
        uptake = fluxes[i]
        fluxes[i] = uptake - cellulase/glc_eq # 1 mol of cellulase produces one mol of glucose equivalents
        i +=1
        
    fluxes =np.append(fluxes, cellulase/glc_eq)

    fluxes *= conc_dict["Growth"]

    if dynamic_system_oligo.pbar is not None:
        dynamic_system_oligo.pbar.update(1)
        dynamic_system_oligo.pbar.set_description('t = {:.3f}'.format(t))
    
    
    return fluxes

dynamic_system_oligo.pbar = None

    
def optimize_parameters_inner_problem_oligo(combination,model,media,rxns,y0,objective_dir,xylog_glc_eq_dict,OD,alternative_solution=False):
    print(combination)
    ts = np.linspace(OD.iloc[0,0], 150, 1000)  
    #vmax_inner_option, km_inner_option, vmax_outer_option, km_outer_option,ki_option = combination

    sol = solve_ivp(
        fun=dynamic_system_oligo,
        t_span=(ts.min(), ts.max()),
        y0=y0,
        t_eval=ts,
        method='LSODA',
        args = (model,rxns,objective_dir,xylog_glc_eq_dict,combination)
    )

    interp_func = interp1d(sol.t, sol.y.T[:,0], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(OD.x.values)

    score_spec = mean_squared_error(OD[" y"].values,y_interp)/(OD[" y"].mean()) 
    penalty =1000*(score_spec)

#    penalty =1000*(1 - r2_score_spec) + 10 * (np.std(y_interp)-np.std(OD[" y"]))**2
    print(f"penalty: {penalty}")    
    
    if alternative_solution:
        return sol,penalty
    else:
        return penalty
    
    
    
def dynamic_system_simplified(t, y,model,rxns,objective_dir,glc_eq_dict,combination,simplified=False):
    """Calculate the time derivative of external species."""    
    
    if simplified==1:
        vmax_outer, Km_outer,Ki = combination
        vmax_inner = 0.18183559
        Km_inner = 2.43662732
        
    elif simplified == 2:
        Km_inner,vmax_outer, Km_outer,Ki = combination
        vmax_inner = 0.18183559
        
    elif simplified ==3:
        Km_inner,vmax_outer, Km_outer = combination
        vmax_inner = 0.18183559
        
        Ki=0
        
    else:
        vmax_inner, Km_inner, vmax_outer, Km_outer,Ki = combination

    rxns_map = copy.copy(rxns)
    rxns_map.remove("EX_xylogluc_e")
    conc_dict = dict(zip(rxns,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        cellulase = add_dynamic_bounds(model, conc_dict,glc_eq_dict,vmax_inner,Km_inner,vmax_outer,Km_outer,Ki)
        
        feasibility = cobra.util.fix_objective_as_constraint(model)
        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns_map, objective_dir)
        
    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    fluxes =lex_constraints.values
    i = 1 
    for key, glc_eq in glc_eq_dict.items():
        uptake = fluxes[i]
        fluxes[i] = uptake - cellulase/glc_eq # 1 mol of cellulase produces one mol of glucose equivalents
        i +=1
        
    fluxes =np.append(fluxes, cellulase)
    fluxes *= conc_dict["Growth"]
    
    return fluxes


def optimize_parameters_inner_problem_simplified(combination,model,media,rxns,y0,objective_dir,xylog_glc_eq_dict,OD,alternative_solution=False,simplified=False):
    print(combination)
    ts = np.linspace(OD.iloc[0,0], 70, 1000)  
    #vmax_inner_option, km_inner_option, vmax_outer_option, km_outer_option,ki_option = combination

    sol = solve_ivp(
        fun=dynamic_system_simplified,
        t_span=(ts.min(), ts.max()),
        y0=y0,
        t_eval=ts,
        method='LSODA',
        args = (model,rxns,objective_dir,xylog_glc_eq_dict,combination,simplified)
    )

    interp_func = interp1d(sol.t, sol.y.T[:,0], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(OD.x.values)

    
    score_spec = mean_squared_error(OD[" y"].values,y_interp)/(OD[" y"].mean()) 
    penalty =1000*(score_spec)

#    penalty =1000*(1 - r2_score_spec) + 10 * (np.std(y_interp)-np.std(OD[" y"]))**2
    print(f"penalty: {penalty}")
      
    
    if alternative_solution:
        return sol,penalty
    else:
        return penalty
    
    
def optimize_with_powell(starting_point,arguments):
    print(f"Started optimizing for {starting_point}")
    result = minimize(optimize_parameters_inner_problem_oligo,
                      starting_point,
                      method="Powell",
                      bounds=[(0.1,1),(0.1,10)],
                      args=arguments,
                     options ={'maxiter':100,'xtol': 1e-4,'ftol':1e-4})
    
    return {
        'starting_point': starting_point,
        'solution': result.x,
        'objective_value': result.fun}
