import sys
sys.path.append("../functions/")

from kinetic_law import michaelis_menten_kinetics, KineticLaw


def cellulase_kinetics(C,met,vmax=2.9,Km=4.4,Ki=11):
    C_spec = C[met]
    bound = -vmax*(C_spec[-1]/((1 + (C["R_EX_cellb_e"][-1]/Ki))*Km + C_spec[-1]))
    return bound


def cellobiose_kinetics(C,met,vmax=5.01,Km=20.22):
    C_spec = C[met]
    bound_mm = michaelis_menten_kinetics(C,"R_EX_cellb_e",vmax=vmax,Km=Km)
    
    
    bound =  max(bound_mm,-0.76)
    #print(f"Bound for cellobiose: {bound}")
    
    return bound



def cellobiose_kinetics_special(C,met,vmax=0.76,Km=0.1,cellulase_kinetic_law_obj=KineticLaw("R_EX_cellulose",kinetic_law_func=cellulase_kinetics)):
    C_spec = C[met]
    cellulase = cellulase_kinetic_law_obj.find_bound(C)
    #print(f"Cellobiose bounds -> Cellulase activity: {round(cellulase,4)}")
    bound_cellulase =  max(0.35*cellulase,-0.76)

    if C_spec[-1]>1e-6:
        bound_mm = michaelis_menten_kinetics(C,"R_EX_cellb_e",vmax=0.76,Km=0.1) ### ??? Apparently R. cellulolyticum has poor control of cellbiose uptake =>vmax is "high" and Km is low?
        bound = min(bound_cellulase,bound_mm)
        
        #print(f"1. Cellulose_conc: {C['R_EX_cellulose_e'][-1]}, Cellobiose_conc: {C['R_EX_cellb_e'][-1]} Cellulase. {cellulase}, bound {bound}")
    else:
        bound = bound_cellulase
        #print(f"2. Cellulose_conc: {C['R_EX_cellulose_e'][-1]}, Cellobiose_conc: {C['R_EX_cellb_e'][-1]} Cellulase. {cellulase}, bound {bound}")
    
    return bound

def glucose_kinetics_special(C,met,vmax,Km,cellulase_kinetic_law_obj=KineticLaw("R_EX_cellulose",kinetic_law_func=cellulase_kinetics)):
    C_spec = C[met]
    cellulase = cellulase_kinetic_law_obj.find_bound(C)
    #print(f"Glucose bounds -> Cellulase activity: {round(cellulase,4)}")
    
    bound =  0.3*cellulase
    
    return bound