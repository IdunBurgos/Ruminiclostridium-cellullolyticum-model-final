# Ruminiclostridium-cellullolyticum-model-final

This is the code for the reconstruction, testing and curation of the R. cellulolyticum strain H10.

It is composed of the following parts

- **1.1.1: Growth on simple sugars.ipynb:** (NB: This is essentially the same as the 'Test of model' notebooks, before it was formalized) In this notebook we get a first glance of how the R. cellulolyticum model is working, in comparison to experimental data.
- **1.1.2. Test for model** (all of these are the same)
- **2.1.1: Curation of model based - cofactors of glycolysis fermentation pathways.ipynb:** In this notebook we focus on incorporating known details of co-factors and metabolic phenotypes in the notebook. 
- **Copy 2.1.1: Curation of model based - cofactors of glycolysis fermentation pathways.ipynb:** A refined version of the Notebook above. (E.g. I added a new reaction HEX1_gtp because the hexokinase reaction was GTP dependent, instead of just changing the stoichiometry (like I did in the Notebook above). 
- **2.1.2: Test-Copy** & **Test**: Essentially the same from what I remember. 
- **2.2.1: Curation of model - False positive gene knockouts.ipynb:** I look into the pathways that are active in the false-positive knockout predictions and try to find the error. 
- **Copy 2.2.1: Curation of model - False positive gene knockouts.ipynb:** Same as above, but like for 'Copy 2.1.1' I tried to be more consistent than I had previously been.
- **2.3.1. Curation of model - Transport and catabolism of galactokinase.ipynb:** Mainly PTS transporters (which also covers galactose metabolism). 
- **3.1.1. Reconstructing pathways for uptake of oligosaccharides from cellulose.ipynb:**
- **3.2.1. Reconstructing pathway for uptake of oligosaccharides from xyloglucan.ipynb**
- **3.3.1. Reconstructing pathways for uptake of oligosaccharides from arabinoxylan.ipynb**
- **Part 4.1. Curation of model based on results from Part 1.1.1 - NGAM and GAM.ipynb**
- **Part 5.1. Removing excessive reactions.ipynb**
