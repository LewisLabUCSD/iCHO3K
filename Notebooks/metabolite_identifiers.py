import pandas as pd
import numpy as np
import json
import urllib.parse as urllib
import requests
import re
from urllib.error import HTTPError


# VARIABLES
ommitSyn = ['MCULE-','SMR','MLS','AKOS','SR-','HMS','EU','OPERA','OPREA','MAYBRIDGE','ZINC','IDI']

def getPubchemCID(cmp,smiles):
    '''
    Retrieves the PubchemID from the Pubchem database

    Input: cmp: Compound Name, smiles
    Output: List of Pubchem IDs
    '''
    query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smiles.replace('/','.'))+'/cids/JSON?MaxRecords=20'
    result = requests.get(query).json()
    if('Fault' in result and cmp!=''):
        query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.quote(cmp)+'/cids/JSON?MaxRecords=20'
        result = requests.get(query).json()
    if (bool(result) and 'IdentifierList' in result):
        ids = list(map(str, filter(lambda a: a != 0, result['IdentifierList']['CID'])))
        return ids
    else:
        return None

def getEBIID(cmp):
    '''
    Retrieves a list of publications regarding a compound of interest

    Input: Compound Name 
    Output: List of publications
    '''
    import re
    full_list = []
    tmp = []
    for item in cmp.split('\n'):
        try:
            val = [s for s in ommitSyn if item[:re.search(r"\d", item).start()].upper() in s]
        except:
            val = []

        if (len(val) == 0):
            query = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&query='+'"'+str(item)+'"'
            result = requests.get(query).json()
            print(result)
            if (bool(result) and ('errCode' not in result) and len(result.keys())>1):
                if(result['hitCount']>0):
                    for i in result['resultList']['result']:
                        if (str(i['id']) not in s for s in full_list): full_list.append((item+'\tPMID'+str(i['id'])+'\t https://pubmed.ncbi.nlm.nih.gov/'+str(i['id']) if str(i['id']).isdigit() else item+'\t'+str(i['id'])))

    return ('\n'.join(full_list))
    

def getChEMBLID(cmp):
    '''
    Retrieves the Chembl ID from the ChEMBL Database

    Input: Compound Name
    Output: ChEMBL ID
    '''
    query = 'https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup/search.json?q='+str(cmp)
    result = requests.get(query).json()
    
    if result and 'chembl_id_lookups' in result:
        chembl_id_lookups = result['chembl_id_lookups']
        if chembl_id_lookups:
            return chembl_id_lookups[0]['chembl_id']
    return None


def getChEMBLSmiID(smi,simil_perc):
    '''
    This function is designed to search for molecules that are similar to a given molecule. 
    It uses the ChEMBL API to perform this search.

    Inputs
    smi : SMILES string representing a molecule.
    simil_perc: Similarity percentage to use for the search (number from 40-100).

    Output
    String that consists of ChEMBL IDs of the similar molecules.

    '''
    ids = []
    query = 'https://www.ebi.ac.uk/chembl/api/data/similarity/'+urllib.quote(smi)+'/'+str(simil_perc)+'.json'
    result = requests.get(query).json()
    #print(result)
    if (bool(result) and 'molecules' in result):
        for mol in result['molecules']:
            if mol['molecule_chembl_id'] not in ids: ids.append(mol['molecule_chembl_id'])
        #print(result['molecule_hierarchy'])
    #return (str(';'.join(ids)))
    return (str('\n'.join(ids)))


def getChEMBLSimilPerc(smi,treshold_simil):
    '''
    Search for molecules that are similar to a given molecule.

    Inputs
    smi : SMILES string representing a molecule.
    treshold_simil: Similarity percentage to use for the search (number from 40-100).

    Output
    List of lists, where each inner list represents a molecule and contains: 
    the ChEMBL ID of the molecule, its canonical SMILES string, and the percentage of similarity to the input molecule.
    '''

    query= 'https://www.ebi.ac.uk/chembl/api/data/similarity/'+urllib.quote(smi)+'/'+str(treshold_simil)+'.json'
    result = requests.get(query).json()
    #print(result)
    if (bool(result) and 'molecules' in result):
        for mol in result['molecules']:

            if (mol['molecule_structures']['canonical_smiles'] != smi and int(float(mol['similarity'])) < 100):
                pairs.append([mol['molecule_chembl_id'],mol['molecule_structures']['canonical_smiles'],mol['similarity']])
    return pairs

def getCIDSmilesInChI(cid):
    '''
    Takes a PubChem ID and returns the Canonical Smiles and InChI for that compound
    '''
    query='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+urllib.quote(str(cid))+'/property/CanonicalSmiles,InChI/JSON'
    result = requests.get(query).json()

    smiles = result['PropertyTable']['Properties'][0]['CanonicalSMILES']
    inchi = result['PropertyTable']['Properties'][0]['InChI']

    return(smiles, inchi)

def getCIDFormula(cid):
    '''
    Takes a PubChem ID and returns the Molecular Formula for that compound
    '''
    query='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+urllib.quote(str(cid))+'/property/MolecularFormula/JSON'
    result = requests.get(query).json()

    formula = result['PropertyTable']['Properties'][0]['MolecularFormula']

    return formula

def getPubChemSubstructure(smi,treshold_simil):
    #[mol ID,canonical smiles]
    substruc = []
    query= 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/'+urllib.quote(smi)+'/cids/JSON?Threshold='+str(treshold_simil)+'&MaxRecords=100'
    result = requests.get(query).json()
    if (bool(result) and 'IdentifierList' in result):
        for p in result['IdentifierList']['CID']:
            substruc.append([p,getCIDSmiles(p)])

    return (substruc)

def getSynonym(cmp,smi):
    
    syn_list = []
    query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smi.replace('/','.'))+'/synonyms/json'
    result = requests.get(query).json()

    # In case SMILES issue but valid, 
    if('Fault' in result and cmp!=''):
        query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.quote(cmp)+'/synonyms/json'
        result = requests.get(query).json()
    if (bool(result) and 'InformationList' in result):
        for item in result['InformationList']['Information']:
            for syn in item['Synonym']:
                if syn not in syn_list and syn != cmp : syn_list.append(syn)

    return str('\n'.join(syn_list))


def getMW(smi):
    '''
    Takes a SMILES string and returns the Molecular Weight for that compound from the PubChem database.
    '''
    query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smi)+'/property/MolecularWeight/json'
    result = requests.get(query).json()
    if (bool(result) and 'PropertyTable' in result):
        result = result['PropertyTable']['Properties'][0]
        return (str(result['MolecularWeight']))
    else:
        return('')


def homogenize_info(df):
    '''
    This function takes the metabolites dataset as input and returns the same dataset with 
    the same information in the "columns_to_homogenize" columns for the same metabolite
    '''
    
    # Define the columns we're interested in
    columns_to_homogenize = ["KEGG","CHEBI","ChEMBLID", "PubChem", "Inchi","EHMNID" ,"SMILES","INCHI2", "CID_ID", "PDB (ligand-expo) Experimental Coordinates  File Url", "Pub Chem Url" ,"ChEBI Url"]
    
    # Create a temporary column for grouping
    df['group_key'] = df['BiGG ID'].str[:-2]

    # Filter out rows with 'not found'
    filtered = df.replace('NaN', np.nan)

    # Group the DataFrame by "BiGG ID" and apply mode function to the groups
    most_frequent = filtered.groupby('group_key')[columns_to_homogenize].agg(lambda x: pd.Series.mode(x.dropna()).iat[0] if not x.dropna().empty else np.nan).reset_index()

    # Merge the most_frequent DataFrame with the original one and update the original DataFrame
    for col in columns_to_homogenize:
        df.set_index('group_key', inplace=True)
        df.update(most_frequent.set_index('group_key')[col])
        df.reset_index(inplace=True)
        
    # Drop the 'group_key' column
    df.drop(columns=['group_key'], inplace=True)

    return df