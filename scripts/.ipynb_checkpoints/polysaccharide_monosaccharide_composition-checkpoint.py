import pandas as pd
import numpy as np
import reframed
from collections import OrderedDict
import molmass

import sys
sys.path.append("../functions/")
import dfba_cobra_multiple_polysaccharides

media = pd.read_csv('../input/media_db.tsv',sep='\t')
media = media.groupby('medium').agg({'compound': set})['compound'].to_dict()

model = dfba_cobra_multiple_polysaccharides.read_model(media)

# The oligosaccharide composition of each polysaccharide

glc_eq_poly_dict =OrderedDict({"EX_arabinoxylan_e":{
"EX_AX_e": 15/6,
"EX_AXX_e": 20/6,
"EX_A23XX_e": 25/6,
"EX_XAXX_e": 25/6,
"EX_XA23XX_e": 30/6},
                              
"EX_xyloglucan_e":{                              
"EX_QLQG_e": 45/6,
"EX_QQLG_e":45/6,
"EX_QLLG_e":51/6,
"EX_QQQG_e":39/6,
"EX_GQQG_e":34/6},

"EX_cellulose_e":{"EX_cellb_e": 2,
"EX_cell3_e": 3,
"EX_cell4_e": 4,
"EX_cell5_e":5},
                               
"EX_xylan_e":{
"EX_xylb_e":10/6,
"EX_xyl3_e":15/6,
"EX_xylan4_e":26/6,
"EX_xylan8_e":52/6}

})

name2id = {"Glucose":"glc__D_e",
                               "Xylose":"xyl__D_e",
                               "Galactose":"gal_e",
                               "Arabinose":"arab__L_e",
                              "Glucuronic acid":"glcur_c"}

oligosac_fracs = {polysac_id:{rxn_id:1/(glc_eq*len(glc_eq_dict)) for rxn_id,glc_eq in glc_eq_dict.items()} for polysac_id,glc_eq_dict in glc_eq_poly_dict.items()}
oligosac_fracs_series = pd.Series({rxn:frac for polysac,dict_ in oligosac_fracs.items() for rxn,frac in dict_.items()})


# The monosaccharide composition of each oligosaccharide
xyloglucan = pd.read_excel("/Users/idunmariaburgos/Documents/polysacc_comp.xlsx",sheet_name="xyloglucan")
cellulose = pd.read_excel("/Users/idunmariaburgos/Documents/polysacc_comp.xlsx",sheet_name="cellulose")
arabinoxylan = pd.read_excel("/Users/idunmariaburgos/Documents/polysacc_comp.xlsx",sheet_name="arabinoxylan")
xylan = pd.read_excel("/Users/idunmariaburgos/Documents/polysacc_comp.xlsx",sheet_name="xylan")

# Process data on monosaccharide composition in each oligosaccharide
monosac_comp = pd.concat([xyloglucan,cellulose,arabinoxylan,xylan])
monosac_comp.fillna(0,inplace=True)

monosac_comp = monosac_comp[monosac_comp.Identifier.str.contains("_e")].copy().drop("Glycosidic bonds",axis=1).set_index("Identifier")


monosac_comp.index = monosac_comp.index.map(lambda x: x.replace("M_","EX_"))
oligosac_fracs_series = oligosac_fracs_series[oligosac_fracs_series.index.isin(monosac_comp.index)]
monosac_comp_mmol = monosac_comp.mul(oligosac_fracs_series,axis=0)

met_molmass = pd.Series({met_name:molmass.Formula(model.metabolites.get_by_id(met_id).formula).mass 
       for met_name,met_id in name2id.items()})

# Convert from mol fraction to weight fraction of monos. in oligos.
monosac_comp_g = monosac_comp_mmol.mul(met_molmass,axis=1)

polysac_molar_weight = pd.DataFrame({polysac_id:monosac_comp_g[monosac_comp_g.index.isin(oligosac_fracs[polysac_id].keys())].sum().to_dict() 
                     for polysac_id,rxns_dict in oligosac_fracs.items()})


polysac_weight_frac = polysac_molar_weight/polysac_molar_weight.sum()

polysac_weight_frac.to_csv("../input/polysaccharide_monosaccharide_weight_fractions.tsv",index=True,sep="\t")