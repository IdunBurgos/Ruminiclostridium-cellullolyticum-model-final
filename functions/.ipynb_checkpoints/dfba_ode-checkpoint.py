from reframed import pFBA, Environment, FVA, FBA, load_cbmodel
import pandas as pd
import numpy as np
from molmass import Formula
import copy

def michaelis_menten_kinetics(C,met,vmax=10,Km=5):
    C_spec = C[met]
    
    bound = -vmax*(C_spec/(C_spec+Km))
    
    return bound

class KineticLaw():
    def __init__(self,met,vmax=10,Km=5,kinetic_law_func=michaelis_menten_kinetics, **kwargs):
        self.met = met
        self.vmax=vmax
        self.Km=Km
        self.kinetic_law_func = kinetic_law_func
        self.extra_params = kwargs
        
    def find_bound(self,C,**kwargs):
        
        all_kwargs = {**self.extra_params, **kwargs}
        
        return self.kinetic_law_func(C=C,met=self.met,vmax=self.vmax,Km=self.Km,**all_kwargs)


    
def cellulase_kinetics(C,met,vmax=2.9,Km=4.4,Ki=11):
    max_tstep = 0.1
    
    C_spec = C[met]
   
    bound = -vmax*(C_spec/((1 + (C["R_EX_cellb_e"]/Ki))*Km + C_spec))
    
    max_bound = -C_spec/(max_tstep*C["R_Growth"])   
    
    return max(bound,max_bound)


def cellobiose_kinetics(C,met,vmax=5.01,Km=0.2):
    max_tstep = 0.1
    
    C_spec = C[met]
    bound = max(-vmax*(C_spec/(Km + C_spec)),-0.76)
    max_bound = -C[met]/(max_tstep*C["R_Growth"])
    
    return max(bound,max_bound)

def glucose_kinetics(C,met,vmax=6.01,Km=0.2):
    max_tstep = 0.1
    
    C_spec = C[met]
    bound = -vmax*(C_spec/(Km + C_spec))
    max_bound = -C[met]/(max_tstep*C["R_Growth"])             
    
    return max(bound,max_bound)


def dynamic_bounds(C,kinetic_laws):
    
    bounds = {}
    
    for met in kinetic_laws.keys():
        bounds[met]=kinetic_laws[met].find_bound(C)
    
    return bounds,kinetic_laws



def cellulase_flux(fluxes,bounds,rxn):
    cellulase =bounds[rxn]
    pfba_cellb = fluxes["R_EX_cellb_e"]
    fluxes["R_EX_cellb_e"] = pfba_cellb - 0.35*cellulase
    pfba_glc = fluxes["R_EX_glc__D_e"]
    fluxes["R_EX_glc__D_e"] = pfba_glc - 0.3*cellulase
    fluxes[rxn]=cellulase
    
    return fluxes



def flux_predictions(t,y,mets,model,kinetic_laws,env,biomass_id,external_rxns_fluxes):

    """
    y = numpy array with concentrations
    mets = array with strings for mets matching the entries in y
    """
    y = np.clip(y, 0, None)
    C = dict(zip(mets,y))
    
    bounds,kinetic_laws = dynamic_bounds(C=C,kinetic_laws=kinetic_laws)
    
    for met,value in bounds.items():
        env[met]=(value,1000)
    env.apply(model,inplace=True,exclusive=True,warning=False)
    
    
    sol = FBA(model)
    
    if sol.status.value == "Optimal":
        fluxes = {rxn:flux for rxn,flux in sol.values.items() if rxn in mets}
        
    else:
        fluxes = {rxn:0 for rxn in model.reactions.keys()}
    
    
    for rxn,func in external_rxns_fluxes.items():
        # Update fluxes based on additional 
        fluxes = func(fluxes,bounds,rxn)
    fluxes_values = np.array([fluxes[rxn] for rxn in C.keys()])*C[biomass_id]
    
    
    
    if flux_predictions.pbar is not None:
            flux_predictions.pbar.update(1)
            flux_predictions.pbar.set_description('t = {:.3f}'.format(t))

    return fluxes_values


    
    
    
    
    
    
    
    