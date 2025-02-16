import cobra

model_cobra = cobra.io.read_sbml_model('../models/RcH10_v2.xml')

check = cobra.manipulation.check_mass_balance(model_cobra)

print(f"Number of imbalanced reactions before: {len(check)}")

model_cobra.metabolites.get_by_id("strcoa_c").charge = -4
model_cobra.metabolites.get_by_id("2agpg161_c").charge = -1
model_cobra.metabolites.get_by_id("2agpg161_p").charge = -1
model_cobra.metabolites.get_by_id("ACP_c").charge = 0
model_cobra.metabolites.get_by_id("pa120_c").charge = -2
model_cobra.metabolites.get_by_id("pa141_c").charge = -2
model_cobra.metabolites.get_by_id("tdeACP_c").charge = 0
model_cobra.metabolites.get_by_id("1hdec9eg3p_c").charge = -2
model_cobra.metabolites.get_by_id("pa161_c").charge = -2
model_cobra.metabolites.get_by_id("ddcaACP_c").charge = 0
model_cobra.metabolites.get_by_id("feenter_c").charge = 3
model_cobra.metabolites.get_by_id("2ahbut_c").charge = -1
model_cobra.metabolites.get_by_id("oc2coa_c").charge = -4
model_cobra.metabolites.get_by_id("pppi_c").charge = -4
model_cobra.metabolites.get_by_id("tsul_c").charge = -2
model_cobra.metabolites.get_by_id("td2coa_c").charge = -4
model_cobra.metabolites.get_by_id("appl_c").charge = 1
model_cobra.metabolites.get_by_id("ametam_c").charge = 2
model_cobra.metabolites.get_by_id("al26da_c").charge = -2
model_cobra.metabolites.get_by_id("r5p_c").charge = -2
model_cobra.metabolites.get_by_id("dcamp_c").charge = -4
model_cobra.metabolites.get_by_id("2shchc_c").charge = -2
model_cobra.metabolites.get_by_id("sbzcoa_c").charge = -5
model_cobra.metabolites.get_by_id("3optcoa_c").charge = -4
model_cobra.metabolites.get_by_id("dscl_c").charge = -7
model_cobra.metabolites.get_by_id("uamr_c").charge = -3
model_cobra.metabolites.get_by_id("udcpdp_e").charge = -3
model_cobra.metabolites.get_by_id("udcpdp_c").charge = -3
model_cobra.metabolites.get_by_id("udcpp_e").charge = -2
model_cobra.metabolites.get_by_id("udcpp_c").charge = -2
model_cobra.metabolites.get_by_id("udcpdp_c").charge = -3
model_cobra.metabolites.get_by_id("g3p_c").charge = -2
model_cobra.metabolites.get_by_id("gtca3_45_BS_c").charge = -45
model_cobra.metabolites.get_by_id("gtca2_45_BS_c").charge = -45
model_cobra.metabolites.get_by_id("gtca1_45_BS_c").charge = -45
model_cobra.metabolites.get_by_id("istfrnA_e").charge = 2
model_cobra.metabolites.get_by_id("2dr1p_c").charge = -2
model_cobra.metabolites.get_by_id("prbamp_c").charge = -4
model_cobra.metabolites.get_by_id("tdecoa_c").charge = -4
model_cobra.metabolites.get_by_id("prbatp_c").charge = -6
model_cobra.metabolites.get_by_id("pphn_c").charge = -2
#model_cobra.metabolites.get_by_id("3hasp__L_c").charge = -1
model_cobra.metabolites.get_by_id("serglugly_c").charge = -1
model_cobra.metabolites.get_by_id("serglugly_e").charge = -1
model_cobra.metabolites.get_by_id("dhpmp_c").charge = -2
model_cobra.metabolites.get_by_id("cspmd_c").charge = 2
model_cobra.metabolites.get_by_id("nadhx__S_c").charge = -2
model_cobra.metabolites.get_by_id("nadhx__R_c").charge = -2
model_cobra.metabolites.get_by_id("ppgpp_c").charge = -6
model_cobra.metabolites.get_by_id("5aizc_c").charge = -3
model_cobra.metabolites.get_by_id("argsuc_c").charge = -1
model_cobra.metabolites.get_by_id("2me4p_c").charge = -2
model_cobra.metabolites.get_by_id("air_c").charge = -2
model_cobra.metabolites.get_by_id("acmanap_c").charge = -2
model_cobra.metabolites.get_by_id("murein3px3p_p").charge = -4
model_cobra.metabolites.get_by_id("murein5px4p_p").charge = -4
model_cobra.metabolites.get_by_id("murein4px4p_p").charge = -4
model_cobra.metabolites.get_by_id("murein5px4px4p_p").charge = -6
model_cobra.metabolites.get_by_id("anhgm3p_p").charge = -2
model_cobra.metabolites.get_by_id("uaagmda_c").charge = -4
model_cobra.metabolites.get_by_id("3hocoa_c").charge = -4

check = cobra.manipulation.check_mass_balance(model_cobra)
print(f"Number of imbalanced reactions now: {len(check)}")

cobra.io.write_sbml_model(model_cobra,'../models/RcH10_v2.xml')