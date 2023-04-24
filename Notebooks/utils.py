from Bio import Entrez
import pandas as pd

def get_gene_info(gene_id):
    '''
    This function retrieves information from the NIH Gene database
    
    input: Gene Entrez ID
    output: Gene Symbol, Gene Name, Gene Description, Gene Ensembl ID, NCBI Transcript ID, NCBI Protein ID
    '''
    Entrez.email = 'account1@theta-ocean-377718.iam.gserviceaccount.com'
    handle = Entrez.efetch(db='gene', id=gene_id, retmode='xml')
    record = Entrez.read(handle)[0]

    gene_name = record['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
    gene_symbol = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']

    # check for different possible formats of the data
    if 'Entrezgene_comments' in record and 'Gene-commentary_comment' in record['Entrezgene_comments'][0]:
        cho_gene_description = record['Entrezgene_comments'][0]['Gene-commentary_comment'][0]['String']
    elif 'Entrezgene_summary' in record:
        cho_gene_description = record['Entrezgene_summary']
    else:
        cho_gene_description = None

    if 'Entrezgene_track-info' in record:
        picr_ensembl_id = next((xref['Dbtag_tag']['Object-id']['Object-id_str'] for xref in record['Entrezgene_gene']['Gene-ref']['Gene-ref_db'] if xref['Dbtag_db'] == 'Ensembl'), None)
    else:
        picr_ensembl_id = None

    xrefs = record['Entrezgene_locus'][0]['Gene-commentary_products']
    
    for xref in xrefs:
        if xref.get('Gene-commentary_accession').startswith('NM_'):
            mRNA_ncbi_id = xref.get('Gene-commentary_accession')
            protein_ncbi_id = xref['Gene-commentary_products'][0].get('Gene-commentary_accession')
            break
        elif xref.get('Gene-commentary_accession').startswith('XM_'):
            mRNA_ncbi_id = xref.get('Gene-commentary_accession')
            protein_ncbi_id = xref['Gene-commentary_products'][0].get('Gene-commentary_accession')
            break
            
    df = pd.read_csv('orthologs&GO.txt')
    for i,row in df.iterrows():
        if row['CHO GeneID'] == gene_id:
            go_terms = row['GO_ids']
            chok1gs_ensembl_id = row['CHO Ensembl ID']
            human_ortholog = str(row['Human GeneID']).replace('.0', '')
            
            human_handle = Entrez.efetch(db='gene', id=human_ortholog, retmode='xml')
            human_record = Entrez.read(human_handle)[0]
            if 'Entrezgene_comments' in human_record and 'Gene-commentary_comment' in human_record['Entrezgene_comments'][0]:
                human_gene_description = human_record['Entrezgene_comments'][0]['Gene-commentary_comment'][0]['String']
            elif 'Entrezgene_summary' in human_record:
                human_gene_description = human_record['Entrezgene_summary']
            else:
                human_gene_description = None
            break
            
        else:
            go_terms = None
            chok1gs_ensembl_id = None
            human_gene_description = None
    
    # Gene description
    if cho_gene_description != None and human_gene_description == None:
        gene_description = cho_gene_description
    elif cho_gene_description == None and human_gene_description != None:
        gene_description = human_gene_description
    else:
        gene_description = None
    
    handle.close()

    return gene_symbol, gene_name, gene_description, picr_ensembl_id, chok1gs_ensembl_id, mRNA_ncbi_id, protein_ncbi_id, go_terms


##### ----- Identification of duplicated reactions ----- #####

import cobra
import numpy as np

def duplicated_reactions(model):
    '''
    This functions retrieves a list of tuples containing
    duplicated reactions from the input cobrapy model.

    input: cobrapy model containing reactions and metabolites.
    output: tuples list with duplicated reactions.
    '''

    #Create stoichiomtric matrix
    S = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')
    #Correlation matrix from stoichiomtric matrix
    cor = np.corrcoef(S.values.T)

    #Extract upper triangle from the correlation matrix
    U = np.triu(cor, k=1)

    #List 's' of identical reactions
    s = np.argwhere(U == 1)

    return s

##### ----- FVA ----- #####
# Single Core
def runMinMax_Single(model, end_rxn_index=None):

    '''
    This function runs a single reaction flux variability analysis (FVA) using linear programming optimization.
    The purpose of this analysis is to determine the maximum and minimum flux values 
    that each reaction in the model can achieve while still satisfying mass balance constraints.
    

    Parameters:
    ---------------
   input: metabolic model represented as a COBRApy Model object

   optional: end_rxn_index argument specifying the index of the last reaction to analyze 
   (default is None, meaning all reactions will be analyzed).

    Returns:
    ---------------
    minmax: array
    
    minmax array with the minimum and maximum values for each reaction side-by-side. 

    '''

    # Import necessary libraries
    import numpy as np
    from cobra.io import load_model

    # If an "end_rxn_index" is not specified, then it is set to "num_rxns", meaning that all reactions in the model will be analyzed.
    num_rxns = len(model.reactions)
    start_rxn_index = 0
    if end_rxn_index is None:
        end_rxn_index = num_rxns

    # Initialize counters and arrays for storing the maximum and minimum flux values for each reaction.
    count = 0
    max_vals = np.zeros((end_rxn_index - start_rxn_index, 1))
    min_vals = np.zeros((end_rxn_index - start_rxn_index, 1))

    # setting each reaction as the objective of the model the loop solves for the maximum and minimum flux values using linear programming.
    for i in range(start_rxn_index, end_rxn_index):
        model.objective = {model.reactions[i]: 1}
        sol = model.optimize('maximize')
        max_vals[count, 0] = sol.fluxes[i]
        sol = model.optimize('minimize')
        min_vals[count, 0] = sol.fluxes[i]
        count += 1
        if count % 500 == 0:
            print(f'{count} reactions analized')
    minmax = np.concatenate((min_vals, max_vals), axis=1)
    return minmax

# Multple Core
def runMinMax_multi(model, end_rxn_index=None, num_processes=None):
    import numpy as np
    from cobra.io import load_model
    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm
    num_rxns = len(model.reactions)
    start_rxn_index = 0
    if end_rxn_index is None:
        end_rxn_index = num_rxns
    if num_processes is None:
        num_processes = 1
    
    reaction_objs = [{} for i in range(start_rxn_index, end_rxn_index)]
    for i in range(start_rxn_index, end_rxn_index):
        reaction_objs[i - start_rxn_index] = {model.reactions[i]: 1}

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        max_fluxes = []
        min_fluxes = []
        tasks = []

        for i in tqdm(range(start_rxn_index, end_rxn_index)):
            tasks.append(executor.submit(model.optimize, 'maximize', reaction_objs[i - start_rxn_index]))
            tasks.append(executor.submit(model.optimize, 'minimize', reaction_objs[i - start_rxn_index]))

        for task in tqdm(tasks, total=len(tasks), desc='Calculating fluxes'):
            flux = task.result().fluxes
            if len(max_fluxes) == 0:
                max_fluxes = flux
                min_fluxes = flux
            else:
                max_fluxes = np.vstack((max_fluxes, flux))
                min_fluxes = np.vstack((min_fluxes, flux))

    minmax = np.concatenate((min_fluxes, max_fluxes), axis=1)
    return minmax

##### ----- Dead-End metabolite Detection ----- #####
# It is the same function as in COBRA matlab toolbox
def detect_dead_ends(model):
    import numpy as np
    from cobra.util import create_stoichiometric_matrix
    
    # Set exchange reaction bounds to allow uptake and secretion of all metabolites
    for rxn_exchange in model.exchanges:
        rxn_exchange.bounds = (-1000, 1000)
    
    # Create the stoichiometric matrix using the dense array type
    S = create_stoichiometric_matrix(model, array_type='dense')
    
    # Get the lower and upper bounds of each reaction
    lb = [i.lower_bound for i in model.reactions]
    ub = [i.upper_bound for i in model.reactions]
    
    # Identify reactions that are irreversible (have a lower bound < 0 and an upper bound > 0)
    ltz = [x < 0 for x in lb]
    gtz = [x > 0 for x in ub]
    
    # Remove columns corresponding to irreversible reactions from the stoichiometric matrix
    S = np.hstack((S[:, gtz], -S[:, ltz]))
    
    # Calculate the absolute sum of each row in the stoichiometric matrix
    abssum = np.sum(np.abs(S), axis=1)
    
    # Calculate the absolute sum of all reactants and products in each reaction
    sumabs = np.abs(np.sum(S, axis=1))
    
    # Identify metabolites that participate only as reactants or only as products in reactions
    onlyConsOrProd = sumabs == abssum
    
    # Identify metabolites that participate in only one reaction
    SPres = S != 0
    onlyOneReac = np.sum(SPres, axis=1) == 1
    
    # Identify metabolites that are not exchanged (i.e., not involved in any exchange reactions)
    ExchangedMets = np.zeros((S.shape[0], 1), dtype=bool)
    
    # Identify metabolites that are dead ends (i.e., participate only as reactants or only as products in reactions and are not exchanged)
    dead_met = ((onlyConsOrProd | onlyOneReac) & ~ExchangedMets)[1]
    
    # Return the list of dead end metabolites
    return dead_met



