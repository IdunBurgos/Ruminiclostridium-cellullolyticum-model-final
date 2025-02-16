import reframed
import copy


model_uni = reframed.load_cbmodel("../models/universe_grampos.xml")
model = reframed.load_cbmodel("../models/RcH10_v2.xml")

def build_mini_model(unique_cofactors=True):
    
    mini_model = reframed.CBModel("mini_model")
    
    # Compartments
    mini_model.add_compartment(reframed.Compartment("C_c","cytosol"))
    mini_model.add_compartment(reframed.Compartment('C_e', "extracellular space"))
    mini_model.add_compartment(reframed.Compartment('C_p', "periplasm"))

    
    # General metabolites (cofactors..)
    mini_model.add_metabolite(copy.copy(model.metabolites["M_h_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_co2_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_h2o_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_h2_e"]))

    mini_model.add_metabolite(copy.copy(model.metabolites["M_h_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_co2_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_h2o_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_h2_c"]))

    mini_model.add_metabolite(copy.copy(model.metabolites["M_nad_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_nadh_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_pi_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_atp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_adp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_gtp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_gdp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_coa_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_fdxo_2_2_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_fdxrd_c"]))

    
    # Tranpsport and exchange of general metabolites
    mini_model.add_reaction(copy.copy(model.reactions["R_H2td"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_h2_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_h_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_CO2t"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_co2_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_H2Ot"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_h2o_e"]))
    
    
    #Specific metabolites
    mini_model.add_metabolite(copy.copy(model.metabolites["M_lac__L_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_etoh_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_etoh_p"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_ac_e"]))
    
    mini_model.add_metabolite(copy.copy(model.metabolites["M_glc__D_e"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_glc__D_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_g6p_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_f6p_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_fdp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_dhap_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_g3p_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_13dpg_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_3pg_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_2pg_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_pep_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_pyr_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_lac__L_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_accoa_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_acald_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_etoh_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_actp_c"]))
    mini_model.add_metabolite(copy.copy(model.metabolites["M_ac_c"]))
    
    #Transport of specific metabolites
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_glc__D_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_ac_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_etoh_e"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_EX_lac__L_e"]))
    
    mini_model.add_reaction(copy.copy(model.reactions["R_ETOHtrpp"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ETOHtex"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_GLCabc"]))
    
    # Glucose to pyruvate
    if unique_cofactors:
        mini_model.add_reaction(copy.copy(model.reactions["R_HEX1_gtp"]))
    else:
        mini_model.add_reaction(copy.copy(model_uni.reactions["R_HEX1"]))
    
    
    if unique_cofactors: 
        mini_model.add_reaction(copy.copy(model.reactions["R_PFK_ppi"]))
    else:
        mini_model.add_reaction(copy.copy(model_uni.reactions["R_PFK"]))
        
    mini_model.add_reaction(copy.copy(model.reactions["R_PGI"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_FBA"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_TPI"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_GAPD"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_PGK"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_PGM"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ENO"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_PYK"]))
    
    
    # Lactate
    mini_model.add_reaction(copy.copy(model.reactions["R_LDH_L"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_L_LACt3"]))
    
    # Ethanol
    mini_model.add_reaction(copy.copy(model.reactions["R_POR_syn"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ACALD"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ALCD2x"]))

    #Acetate
    mini_model.add_reaction(copy.copy(model.reactions["R_PTAr"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ACKr"]))
    #mini_model.add_reaction(copy.copy(model.reactions["R_ACt2r_1"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ACt2r"]))
    
    # Cofactor balancing
    mini_model.add_reaction(copy.copy(model.reactions["R_FNRR"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_FNRR3"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ATPM"]))
    mini_model.add_reaction(copy.copy(model.reactions["R_ATPS4r"]))
    

    mini_model.add_reaction_from_str("R_GTPM: M_gtp_c + M_h2o_c --> M_gdp_c + M_h_c + M_pi_c")
    return mini_model