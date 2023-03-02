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

