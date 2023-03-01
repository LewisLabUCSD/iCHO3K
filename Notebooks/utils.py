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
        gene_description = record['Entrezgene_comments'][0]['Gene-commentary_comment'][0]['String']
    elif 'Entrezgene_summary' in record:
        gene_description = record['Entrezgene_summary']
    else:
        gene_description = None

    if 'Entrezgene_track-info' in record:
        gene_ensembl_id = next((xref['Dbtag_tag']['Object-id']['Object-id_str'] for xref in record['Entrezgene_gene']['Gene-ref']['Gene-ref_db'] if xref['Dbtag_db'] == 'Ensembl'), None)
    else:
        gene_ensembl_id = None

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
    

    handle.close()

    return gene_symbol, gene_name, gene_description, gene_ensembl_id, mRNA_ncbi_id, protein_ncbi_id