import os, re;
import gzip, time;
import pandas as pd;
import multiprocessing;

outdir="/path/to/workdir"
os.chdir(outdir)

yeastData=pd.read_csv('https://raw.githubusercontent.com/cpdong/imParalog/main/data/yeast_essential_hsapiens_orthologs.tsv',header=0, sep='\t')
yeastData = yeastData[yeastData["Source"] == 'OGEE v3']
yeast_essential_gene = yeastData['hsapiens_homolog_ensembl_gene_id'].unique()


# Running CodeML module function from PAML
paralog_bioMart_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
para_biomart = pd.read_csv(paralog_bioMart_file, header=0, sep="\t")
para_biomart = para_biomart[["ensembl_gene_id","hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

yeast_essentials = []; 
for index, row in para_biomart.iterrows():
    if any(x in yeast_essential_gene for x in row):
        yeast_essentials.append(1)
    else:
        yeast_essentials.append(0)
        
para_biomart['yeast_essentials'] = yeast_essentials;

para_biomart.to_csv('Paralog_gene_yeast_essentials.tsv', index=False, sep='\t');
#
