import os, re;
import gzip, time;
import subprocess,shutil;
from Bio import SeqIO;
import pandas as pd;
import urllib.request;

os.chdir('/path/to/directory')

# Generate pure gene-transcript-protein-pepseq-cdsseq info
paralog_bioMart_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
paralogData= pd.read_csv(paralog_bioMart_file, header=0, sep="\t")
pdata1=paralogData[["ensembl_gene_id", "hsapiens_paralog_canonical_transcript_protein"]].drop_duplicates(keep='first')
pdata1.columns = ['gene_id', 'protein_id']
pdata2=paralogData[["hsapiens_paralog_ensembl_gene","hsapiens_paralog_ensembl_peptide"]].drop_duplicates(keep='first')
pdata2.columns = ['gene_id', 'protein_id']
paras = pd.concat([pdata1, pdata2], ignore_index = True).drop_duplicates(keep='first')
paras_dict=pd.Series(paras.gene_id.values,index=paras.protein_id).to_dict()

# Peptide data
urllib.request.urlretrieve("https://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", "pcgene_translation.fa.gz")
human_pep_fa="pcgene_translation.fa.gz"
pepData=[['geneid','txid','proteinid','pepSeq']]
open_method = gzip.open if human_pep_fa.lower().endswith('.gz') else open;
with open_method(human_pep_fa, "rt") as fa:
    for s in SeqIO.parse(fa, "fasta"):
        id = s.id[:15]
        if id in paras_dict:
            geneid=[ x.split(':')[1] for x in s.description.split(' ') if 'ENSG' in x ][0][:15]
            txid=[ x.split(':')[1] for x in s.description.split(' ') if 'ENST' in x ][0][:15]
            pepData.append([geneid,txid,id,str(s.seq)])
pep_df = pd.DataFrame(pepData[1:],columns=pepData[0]);

# CDS data
urllib.request.urlretrieve("https://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz", "pcgene_cds.fa.gz")
human_cds_fa="Homo_sapiens.GRCh38.cds.all.fa"
cdsData=[['txid','cdsSeq']]
open_method = gzip.open if human_cds_fa.lower().endswith('.gz') else open;
with open_method(human_cds_fa, "rt") as fa:
    for s in SeqIO.parse(fa, "fasta"):
        txid = s.id[:15]
        geneid=[ x.split(':')[1] for x in s.description.split(' ') if 'ENSG' in x ][0][:15]
        if geneid in paras_dict.values():
            cdsData.append([txid, str(s.seq)]) 
cds_df = pd.DataFrame(cdsData[1:],columns=cdsData[0]);      
pep_cds_df = pd.merge(pep_df,cds_df, on='txid', how="inner")
pep_cds_df.shape
pep_cds_df.to_csv('./Ensembl102_gene_pep_cds_sequence_table.txt', index=False, sep='\t');
