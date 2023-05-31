import json
import urllib.parse as urllib
import requests
import re
from urllib.error import HTTPError


# VARIABLES
ommitSyn = ['MCULE-','SMR','MLS','AKOS','SR-','HMS','EU','OPERA','OPREA','MAYBRIDGE','ZINC','IDI']

def getPubchemCID(cmp,smiles):
    query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smiles.replace('/','.'))+'/cids/JSON?MaxRecords=20'
    result = requests.get(query).json()
    if('Fault' in result and cmp!=''):
        query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.quote(cmp)+'/cids/JSON?MaxRecords=20'
        result = requests.get(query).json()
    if (bool(result) and 'IdentifierList' in result):
        ids = filter(lambda a: a != 0, result['IdentifierList']['CID'])
        return (';'.join(str(x).replace("'",'') for x in ids))
    else:
        return (' ')

def getEBIID(id):
    import re
    full_list = []
    tmp = []
    for item in id.split('\n'):
        try:
            val = [s for s in ommitSyn if item[:re.search(r"\d", item).start()].upper() in s]
        except:
            val = []

        if (len(val) == 0):
            query = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&query='+'"'+str(item)+'"'
            result = requests.get(query).json()
            #print(result)
            if (bool(result) and ('errCode' not in result) and len(result.keys())>1):
                if(result['hitCount']>0):
                    for i in result['resultList']['result']:
                        if (str(i['id']) not in s for s in full_list): full_list.append((item+'\tPMID'+str(i['id'])+'\t https://pubmed.ncbi.nlm.nih.gov/'+str(i['id']) if str(i['id']).isdigit() else item+'\t'+str(i['id'])))

    return ('\n'.join(full_list))
    

def getChemblID(id):
    full_list = []
    query = 'https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup/search.json?q='


def getChemblSmiID(smi,simil_perc):
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


def getChemblSimilPerc(smi,treshold_simil):
    #[mol ID,canonical smiles, similarity percentage]
    pairs = []

    query= 'https://www.ebi.ac.uk/chembl/api/data/similarity/'+urllib.quote(smi)+'/'+str(treshold_simil)+'.json'
    result = requests.get(query).json()
    #print(result)
    if (bool(result) and 'molecules' in result):
        for mol in result['molecules']:

            if (mol['molecule_structures']['canonical_smiles'] != smi and int(float(mol['similarity'])) < 100):
                pairs.append([mol['molecule_chembl_id'],mol['molecule_structures']['canonical_smiles'],mol['similarity']])
    return pairs

def getCIDSmiles(cid):
    query='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+urllib.quote(str(cid))+'/property/CanonicalSmiles/JSON'
    result = requests.get(query).json()
    #print(result) UNCOMMENTED
    id = result['PropertyTable']['Properties'][0]['CanonicalSMILES']

    return(id)


def getPubChemSubstructure(smi,treshold_simil):
    #[mol ID,canonical smiles]
    substruc = []
    #query= 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/'+urllib.quote(smi)+'/cids/JSON?Threshold='+str(treshold_simil)+'&MaxRecords=100'
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
    query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smi)+'/property/MolecularWeight/json'
    #print(query)
    result = requests.get(query).json()
    if (bool(result) and 'PropertyTable' in result):
        result = result['PropertyTable']['Properties'][0]
        #print(result)
        return (str(result['MolecularWeight']))
        #return (str(result['PropertyTable']['Properties']['MolecularWeight']))
    else:
        return('')