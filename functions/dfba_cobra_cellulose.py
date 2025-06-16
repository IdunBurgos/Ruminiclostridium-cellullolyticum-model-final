import copy
import cobra

import numpy as np
import pandas as pd

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error 


cellulose_exp = pd.read_csv("../input/Desvaux2001_batch_data/cellulose_g.csv")
biomass_exp = pd.read_csv("../input/Desvaux2001_batch_data/biomass_mg.csv")

cellulose_exp = pd.read_csv("../input/Desvaux2001_batch_data/cellulose_g.csv")
molar_mass = 173.85  # Based on glucose equivalent (based on 0.35 cellobiose and 0.3 glucose)
cellulose_exp["y mmol"]= cellulose_exp[" y"].map(lambda x: (x/molar_mass)*1000)
cellulose_exp = cellulose_exp[cellulose_exp.x<70].copy()

biomass_exp = pd.read_csv("../input/Desvaux2001_batch_data/biomass_mg.csv")
biomass_exp["y g"] = biomass_exp[" y"].map(lambda x: x/1000)
biomass_exp = biomass_exp[biomass_exp.x<70].copy()

def add_dynamic_bounds(model, conc_dict,combination=[6.01,0.2,5.01,0.2,2.9]):
    """Use external concentrations to bound the uptake flux of glucose."""
    
    vmax_inner_glc, Km_inner_glc,vmax_inner_cellb, Km_inner_cellb,vmax_outer= combination
                                                     
    Km_outer,Ki = [4.4,11]
    
    glucose_max_import = -vmax_inner_glc* conc_dict["EX_glc__D_e"] / (Km_inner_glc + conc_dict["EX_glc__D_e"])
    cellobiose_max_import = -vmax_inner_cellb * conc_dict["EX_cellb_e"]/ (Km_inner_cellb + conc_dict["EX_cellb_e"]) 
    
    cellulase = -vmax_outer*(conc_dict["EX_cellulose_e"]/((1 + (conc_dict["EX_cellb_e"]/Ki))*Km_outer + conc_dict["EX_cellulose_e"]))    

    model.reactions.EX_glc__D_e.lower_bound = glucose_max_import
    model.reactions.EX_cellb_e.lower_bound = cellobiose_max_import
    
    return cellulase
    
    
    
    
def dynamic_system(t, y,model,rxns,objective_dir,combination):
    """Calculate the time derivative of external species."""
    rxns_map = copy.copy(rxns)
    rxns_map.append("EX_cellulose_e")
    
    conc_dict = dict(zip(rxns_map,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        cellulase = add_dynamic_bounds(model, conc_dict,combination)
        
        feasibility = cobra.util.fix_objective_as_constraint(model)

        lex_constraints = cobra.util.add_lexicographic_constraints(model, rxns, objective_dir)

        
    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    fluxes = lex_constraints.values

    glucose_uptake = fluxes[1]
    cellobiose_uptake = fluxes[2]
    
    fluxes[1] = glucose_uptake - cellulase*0.35
    fluxes[2] = cellobiose_uptake - cellulase*0.3
    fluxes =np.append(fluxes, cellulase)
    
    fluxes *= conc_dict["Growth"]

    # This implementation is **not** efficient, so I display the current
    # simulation time using a progress bar.
    if dynamic_system.pbar is not None:
        dynamic_system.pbar.update(1)
        dynamic_system.pbar.set_description('t = {:.3f}'.format(t))
    
    return fluxes

dynamic_system.pbar = None




def optimize_parameters(combination,model,rxns,y0,objective_dir,alternative_solution=False,t_end=False):

    
    try:
        if t_end:
            ts = np.linspace(biomass_exp.iloc[0,0], t_end, 1000)   
        else:
            ts = np.linspace(biomass_exp.iloc[0,0], biomass_exp.iloc[biomass_exp["x"].size-1,0], 1000)   

        sol = solve_ivp(
            fun=dynamic_system,
            t_span=(ts.min(), ts.max()),
            y0=y0,
            t_eval=ts,
            method='LSODA',
            events = [infeasible_event],
            args = (model,rxns,objective_dir,combination)
        )
    except Exception as e:
        print(f"\t had issues with this combination: {combination}\nException: {e}")
        return 1e6 # Returns high penalty score if combination is infeasible...
    rxns_extra = rxns.copy()
    rxns_extra.append("EX_cellulose_e")
    C_dict_results = dict(zip(rxns_extra,sol.y))
    
    # Growth
    
    interp_func = interp1d(sol.t,C_dict_results["Growth"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(biomass_exp.x.values)
    
    score_spec = mean_squared_error(biomass_exp["y g"].values,y_interp)/(biomass_exp["y g"].mean()) 
    penalty_growth =1000*(score_spec)
    
    
    #  Cellulose
    interp_func = interp1d(sol.t, C_dict_results["EX_cellulose_e"], kind='linear', fill_value="extrapolate")
    y_interp = interp_func(cellulose_exp.x)
    
    score_spec = mean_squared_error(cellulose_exp["y mmol"],y_interp)/(cellulose_exp["y mmol"].mean()) 
    penalty_cellulose =10*(score_spec)

    penalty = penalty_growth + penalty_cellulose
    
    print(f"\t penalty: {penalty} for combination: {combination}")
    
    if alternative_solution:
        return sol,penalty,{"penalty_growth":penalty_growth,"penalty_cellulose":penalty_cellulose}
    else:
        return penalty
    


def infeasible_event(t, y,model,rxns,objective_dir,combination):
    """
    Determine solution feasibility.

    Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
    This function re-solves the LP to determine whether or not the solution is feasible
    (and if not, how far it is from feasibility). When the sign of this function changes
    from -epsilon to positive, we know the solution is no longer feasible.
    """
    rxns_map = copy.copy(rxns)
    rxns_map.append("EX_cellulose_e")
    
    conc_dict = dict(zip(rxns_map,y))    
    with model:


        add_dynamic_bounds(model, conc_dict,combination)

        feasibility = cobra.util.fix_objective_as_constraint(model)

    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True