def michaelis_menten_kinetics(C,met,vmax=10,Km=5):
    
    
    C_spec = C[met]
    
    bound = -vmax*(C_spec[-1]/(C_spec[-1]+Km))
    
    if met=="R_EX_glc__D_e":
        print(f"Using michaelis menten kinetics for glucose: {bound}")
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