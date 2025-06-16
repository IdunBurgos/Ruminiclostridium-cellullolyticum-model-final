import pandas as pd
import reframed
import cobra
import re
from reframed import GPRAssociation, Protein

def parse_gpr(gpr_string):
# Made in cooperation with github copilot

    # Clean up the string
    tokens = re.findall(r'\(|\)|\w+|and|or', gpr_string.replace('AND', 'and').replace('OR', 'or'))

    def parse_tokens(tokens):
        stack = []
        current = []
        op = None
        while tokens:
            token = tokens.pop(0)
            if token == '(':
                group = parse_tokens(tokens)
                if op == 'and':
                    if isinstance(current[-1], list):
                        current[-1] += group
                    else:
                        current[-1] = [current[-1]] + group
                elif op == 'or':
                    current += group
                else:
                    current += group
            elif token == ')':
                break
            elif token in ('and', 'or'):
                op = token
            else:
                if op == 'and' and current:
                    if isinstance(current[-1], list):
                        current[-1].append(token)
                    else:
                        current[-1] = [current[-1], token]
                else:
                    current.append(token)
        return current

    def build_gpr_tree(tree):
        # tree can be nested lists (for ands) or flat lists (for ors)
        gpr = GPRAssociation()
        if isinstance(tree, list) and any(isinstance(i, list) for i in tree):
            # OR of ANDs
            for item in tree:
                protein = Protein()
                protein.genes = item if isinstance(item, list) else [item]
                gpr.proteins.append(protein)
        elif isinstance(tree, list):
            # OR of single genes
            for item in tree:
                protein = Protein()
                protein.genes = [item]
                gpr.proteins.append(protein)
        else:
            # Single gene
            protein = Protein()
            protein.genes = [tree]
            gpr.proteins.append(protein)
        return gpr

    parsed = parse_tokens(tokens)
    return build_gpr_tree(parsed)


# Load model
print("Load model '../models/other_models/h10-C_cellulolyticium.xml'")
model_ifs = reframed.load_cbmodel("../models/other_models/h10-C_cellulolyticium.xml",load_gprs=False,flavor="cobra", exchange_detection="^EX_.*")
# Load data

print("Load model data")
ifs_data = pd.read_excel("../input/biot_201000159_sm_suppinfo_2.xls",sheet_name="REACTIONS")
ifs_data.GENE = ifs_data.GENE.fillna("")
ifs_data = ifs_data[ifs_data.GENE.str.match(r"^(Ccel_|\()")]

#Add gpr to iFS431 model
print("Add GPRs to model...")
for rxn, gpr_string in ifs_data.set_index("REACTION ABBREVIATION")["GENE"].iteritems():
    gpr_string = gpr_string.replace("CAC","Ccel_").replace("CAP","Ccel_")
    
    gpr = parse_gpr(gpr_string)
    if not isinstance(rxn,str):
        continue
    
    rxn_id = rxn
    
    if (rxn_id not in model_ifs.reactions): 
    
        continue
        
        
    model_ifs.set_gpr_association(rxn_id,gpr)

model_ifs.remove_metabolites([met for met in model_ifs.metabolites if met.endswith("_b")])
model_ifs.update()

print("Save model as '../models/other_models/iFS431_genes.xml'")
reframed.save_cbmodel(model_ifs,"../models/other_models/iFS431_genes.xml",flavor="fbc")

model_cobra = cobra.io.read_sbml_model("../model/other_models/iFS431_genes2.xml")


