import reframed
import re
from collections import OrderedDict
import copy
import requests
import time

from concurrent.futures import ThreadPoolExecutor, as_completed


## DATA FOR CURATION

map_mets = {"seedm/":"seed.compound:",
        "seed.compound/":"seed.compound:",
       "biggm/":"bigg.metabolite:",
        "bigg.metabolite/":"bigg.metabolite:",
       "metacycm/":"metacyc.compound:",
       "keggc/":"kegg.compound:",
       "kegg.compound/":"kegg.compound:",
       "keggd/":"kegg.drug:",
       "kegg.drug/":"kegg.drug:",
       "chebi/CHEBI:":"CHEBI:",
       "sabiorkm/":"sabiork.compound:",
       "sabiork.compound/":"sabiork.compound:",
       "metacyc.compound/":"metacyc.compound:",
       "lipidmapsm/":"lipidmaps:",
       "lipidmaps/":"lipidmaps:",
       "slm/":"SLM:",
       "reactome/":"reactome:",
       "reactomem/":"reactome:",
       "metanetx.chemical/":"metanetx.chemical:",
       "pubchem.compound/":"pubchem.compound:",
       "hmdb/":"hmdb:",
       "inchikey/":"inchikey:",
       "inchi/":"inchi:"}

map_rxns = {"seedr/":"seed.reaction:",
        "seed.reaction/":"seed.reaction:",
        "biggr/":"bigg.reaction:",
        "bigg.reaction/":"bigg.reaction:",
        "metacycr/":"metacyc.reaction:",
        "keggr/":"kegg.reaction:",
        "kegg.reaction/":"kegg.reaction:",
        "sabiorkr/":"sabiork.reaction:",
        "sabiork.reaction/":"sabiork.reaction:",
        "metacyc.reaction/":"metacyc.reaction:",      
        "reactome/":"reactome:",
        "reactomer/":"reactome:",
        "metanetx.reaction/":"metanetx.reaction:",
       "rhear/":"rhea:",
       "rhea/":"rhea:",
       "brenda/":"brenda:",
       "ec-code/":"ec-code:",
       "biocyc/":"biocyc:"}


M_prefix = ["seed","bigg","kegg"]


def is_resolvable_http(uri, max_retries=5, base_delay=2):

    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36"
    }
    
    timeout = 2  # Initial timeout
    
    for attempt in range(max_retries):
        try:
            response = requests.get(uri, allow_redirects=True, timeout=timeout, headers=headers)
            return uri, response.status_code < 400  # Return tuple (uri, result)
        
        except requests.Timeout:
            if attempt==0:
                print(f"\t Timeout for {uri}. Retrying...")
        
        except requests.ConnectionError:
            if attempt==0:
                print(f"\t Connection error {uri}. Retrying...")

        except requests.RequestException as e:
            print(f"\t Other request error for {uri}: {e}")
            return uri, False  # Return (uri, False)

        # Wait before retrying
        wait_time = base_delay * (2 ** attempt)
        time.sleep(wait_time)

        # Increase timeout but limit it
        timeout = min(timeout * 2, 10)  

    print(f"\t Failed to connect to {uri} after {max_retries} retries.")
    return uri, False  # Return (uri, False)


# Function to run multiple URI checks in parallel

def check_uris_parallel(uris, max_workers=5):
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_uri = {executor.submit(is_resolvable_http, uri): uri for uri in uris}
        
        for future in as_completed(future_to_uri):
            uri, result = future.result()
            results[uri] = result

    return results  # Dictionary {uri: True/False}


def check_uri_validitiy(elements,ismet=True):
    uris = []
    kept_elements = []
    
    for element in elements:
        match = re.search(r'rdf:resource="([^"]+)"',element)
        if match:
            uri = match.group(1)
            uris.append(uri)
            kept_elements.append(element)

    map_uris = dict(zip(uris,kept_elements))
    
    max_workers = min(10,len(uris))
    
    results = check_uris_parallel(uris,max_workers=max_workers)
    
    results_mapped = {map_uris[key]:value for key, value in results.items()}

    
    return results_mapped

def process_metadata(model, met,ismet=True):
    if ismet:
        metadata = model.metabolites[met].metadata
    else:
        metadata = model.reactions[met].metadata
    if "XMLAnnotation" not in metadata:
        return  # Skip if no metadata
    
    metadata_list = metadata["XMLAnnotation"].split("\n")
    elements = [element for element in metadata_list if "<rdf:li rdf:resource=" in element]
    if len(elements)==0:
        return
    
    # Automatically set elements to True
    is_included = OrderedDict((element, True) for element in metadata_list)
    
    # Check all URIs in parallel
    results_mapped = check_uri_validitiy(elements,ismet=ismet)
    
    # Update to new results
    is_included.update(results_mapped)
    
    # Filter metadata based on valid URIs
    metadata_list_new =[uri for uri in is_included.keys() if is_included.get(uri, True)]
    if len(metadata_list_new)<len(metadata_list):
        print(f"\t Removed problematic URIs for {met}")
    
    metadata_new = "\n".join(metadata_list_new)
    return metadata_new


def change_uri(element,ismet=True):
    if ismet:
        map_=map_mets
        prefix = "M_"
    else:
        map_=map_rxns
        prefix = "R_"
    changed = 0

    # Check for and replace elements with the wrong namespace prefix
    for key, value in map_.items():
        if key in element:
            element_new = element.replace(key,value)
            changed+=1
            break

    if changed==0:
        element_new = element        

    # Replace prefix of M_ where necessary for specified database cross-references
    for value in M_prefix:
        if value not in element_new:
            continue
        if f":{prefix}" in element_new:
            element_new = element_new.replace(f":{prefix}",":")
        elif f"/{prefix}" in element_new:
            element_new = element_new.replace(f"/{prefix}",":")
    return element_new



### Fix metabolites

print("loading models...")
model = reframed.load_cbmodel("../models/RcH10_v5.xml")

model_draft = reframed.load_cbmodel("../models/RcH10_draft.xml")

print("fixing metabolites...")
for met in model.metabolites:
    
    metadata = model.metabolites[met].metadata
    if "XMLAnnotation" not in metadata.keys():
        continue
        
    metadata_hasxml = metadata["XMLAnnotation"]
    metadata_list = metadata_hasxml.split("\n")
      
    metadata_list_new = []

    for element in metadata_list:

        # Envipath annotations seem to be incorrect
        if ("envipath" in element) or ("hmdb" in element):
            continue
            
        element_new = change_uri(element,ismet=True)
        
        # split duplicated items in the database
        match = re.search(r'(.*?pubchem\.compound:)([\d\s]+)"', element_new)
        if match:
            prefix = match.group(1)  # Extract the prefix
            numbers = match.group(2).strip().split()  # Extract the numbers, removing extra spaces
            for number in numbers:
                metadata_list_new.append(prefix+number+'"/>' )
        else:
            metadata_list_new.append(element_new)
            

    metadata_list_new = list(OrderedDict.fromkeys(metadata_list_new))
    metadata_new = "\n".join(metadata_list_new)
    model.metabolites[met].metadata["XMLAnnotation"] = metadata_new

model.update()

#### For all relevant URIs: check if URI is valid
print("checking URIs for metabolites...\n Having problems with following URIs:")
for met in model.metabolites:
    metadata_new = process_metadata(model,met)
    model.metabolites[met].metadata["XMLAnnotation"] = metadata_new


print("Fixing reactions...")
for rxn in model.reactions:
    metadata = model.reactions[rxn].metadata
    if "XMLAnnotation" not in metadata.keys():
        continue
    
    # FOR TRANSPORT AND EXCHANGE REACTIONS SET METDATA FROM DRAFT
    reaction_type = model.reactions[rxn].reaction_type
    if (reaction_type==reframed.ReactionType.TRANSPORT) or (reaction_type==reframed.ReactionType.EXCHANGE) :
        
        # If reaction is in draft
        if rxn in model_draft.reactions:
            draft_metadata = model_draft.reactions[rxn].metadata
            
            # SET NEW METADATA
            if "XMLAnnoation" in draft_metadata.keys():
                model.reactions[rxn].metadata["XMLAnnoation"] = copy.copy(model_draft.reactions[rxn].metadata["XMLAnnoation"])   
        
    metadata_hasxml = metadata["XMLAnnotation"]
    metadata_list = metadata_hasxml.split("\n")
    metadata_list_new = []

    for element in metadata_list:
        #Brenda references do not work
        if ("brenda" in element):
            continue
            
        # exchange reactions should not have other associations with other bigg reactions
        if ("bigg" in element) and (rxn[2:].startswith("EX_")) and (rxn[2:] not in element):
            continue
            
        element_new = change_uri(element,ismet=False)
        
        metadata_list_new.append(element_new)
        
    metadata_list_new = list(OrderedDict.fromkeys(metadata_list_new))
    metadata_new = "\n".join(metadata_list_new)
    model.reactions[rxn].metadata["XMLAnnotation"] = metadata_new

model.update()


print("checking URIs for reactions...\n Having problems with following URIs:")
for rxn in model.reactions:
    metadata_new = process_metadata(model,rxn,ismet=False)
    model.reactions[rxn].metadata["XMLAnnotation"] = metadata_new
    
print("save model...")
reframed.save_cbmodel(model,"../models/RcH10_final.xml")