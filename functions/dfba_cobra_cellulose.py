import copy
import cobra

import numpy as np

def add_dynamic_bounds(model, conc_dict):
    """Use external concentrations to bound the uptake flux of glucose."""
    
    glucose_max_import = -6.01* conc_dict["EX_glc__D_e"] / (0.2 + conc_dict["EX_glc__D_e"])
    cellobiose_max_import = max((-5.01 * conc_dict["EX_cellb_e"]/ (0.2 + conc_dict["EX_cellb_e"]) ),-0.76)
    
    cellulase = -2.9*(conc_dict["EX_cellulose_e"]/((1 + (conc_dict["EX_cellb_e"]/11))*4.4 + conc_dict["EX_cellulose_e"]))    

    model.reactions.EX_glc__D_e.lower_bound = glucose_max_import
    model.reactions.EX_cellb_e.lower_bound = cellobiose_max_import
    
    return cellulase
    
    
    
    
def dynamic_system(t, y,model,rxns,objective_dir):
    """Calculate the time derivative of external species."""
    rxns_map = copy.copy(rxns)
    rxns_map.append("EX_cellulose_e")
    
    conc_dict = dict(zip(rxns_map,y))
    
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:
        cellulase = add_dynamic_bounds(model, conc_dict)
        
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


def infeasible_event(t, y,model,rxns,objective_dir):
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


        add_dynamic_bounds(model, conc_dict)

        feasibility = cobra.util.fix_objective_as_constraint(model)

    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True