#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:21:25 2023

@author: cpdong
"""

import os, zipfile;
import pandas as pd
import numpy as np
from pathlib import Path
import urllib.request;

def paralog_closest_check(paralogData_from_biomart):
    # paralogData_from_biomart = paralogData.copy()
    p0 = paralogData_from_biomart[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene', 'hsapiens_paralog_perc_id']]; p0.columns = ['geneid', 'paralog_gene','perc_idt'];
    p1 = paralogData_from_biomart[['ensembl_gene_id', 'hsapiens_paralog_perc_id']]; p1.columns = ['geneid', 'perc_idt'];
    p2 = paralogData_from_biomart[['hsapiens_paralog_ensembl_gene', 'hsapiens_paralog_perc_id']]; p2.columns = ['geneid', 'perc_idt'];
    paralog_max_idt = pd.concat([p1,p2])
    paralog_max_idt = pd.DataFrame(paralog_max_idt.groupby('geneid')['perc_idt'].max()).rename(columns={'perc_idt':'perc_idt_max'})
    
    pdat = pd.merge(p0, paralog_max_idt, on='geneid', how='left')
    pdat['closest'] = np.where(pdat['perc_idt'] == pdat['perc_idt_max'], 1, 0)
    pdat = pdat.drop(columns=['perc_idt','perc_idt_max']);
    
    return pdat
    
def paralog_familysize(paralogPairs):
    #paralogPairs = paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']];
    paralogPairs.columns = ['gene1','gene2']
    paralogPairs_all = pd.concat([paralogPairs,
                                  paralogPairs[['gene2','gene1']].rename(columns={'gene2':'gene1', 'gene1':'gene2'})
                                  ]).drop_duplicates(keep='first')
    
    # Set of all paralogs for each A1
    paralog_families = paralogPairs_all.groupby('gene1').agg({'gene2':set}).reset_index()
    # Family size = union of all A1 and all A2 paralogs
    df = pd.merge(paralogPairs, paralog_families.rename(columns={'gene2':'A1_paralog'}), on='gene1', how='left')
    df = pd.merge(df, paralog_families.rename(columns={'gene1':'gene2', 'gene2':'A2_paralog'}), on='gene2', how='left')
    
    df['A1_paralog'] = df['A1_paralog'].apply(lambda d: d if not pd.isnull(d) else set())
    df['A2_paralog'] = df['A2_paralog'].apply(lambda d: d if not pd.isnull(d) else set())
    
    
    df['family'] = df.apply(lambda x: frozenset(x.A1_paralog.union(x.A2_paralog)), axis=1)
    df['family_size'] = df.apply(lambda x: len(x.A1_paralog.union(x.A2_paralog)), axis=1)
    df = df.drop(columns=['A1_paralog','A2_paralog','family']); # can keep family col if check
    df.columns = ['geneid','paralog_gene','familysize']    
    
    return df

def paralog_chrom_check(paralogPairs):
    #paralogPairs = paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']];
    paralogPairs.columns = ['gene1','gene2']

    annot = pd.read_csv('https://raw.githubusercontent.com/cpdong/ParalogICB/main/data/gencode.v36_annotation.tsv',header=0,sep='\t')
    annot = annot[['Geneid', 'Chrom']]
    
    df = pd.merge(paralogPairs, annot.rename(columns={'Geneid':'gene1', 'Chrom':'A_chrom'}), on='gene1', how='left')
    df = pd.merge(df, annot.rename(columns={'Geneid':'gene2', 'Chrom':'B_chrom'}), on='gene2', how='left')
    
    df['both_same_chrom'] = np.where(df['A_chrom'] == df['B_chrom'], 1, 0)
    df = df.drop(columns=['A_chrom','B_chrom']); # can keep family col if check
    df.columns = ['geneid','paralog_gene','both_same_chrom']
    return df


def paralog_HPA_exprs(paralogPairs, tmpdir='./'):
    if tmpdir=='./':
        tmpdir = os.getcwd()
    # paralogPairs = paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']];
    
    tmpdir='/Users/cpdong/Downloads/'
    hpa_url = 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip'
    urllib.request.urlretrieve(hpa_url, str(Path(tmpdir)) + '/proteinatlas.tsv.zip');
    with zipfile.ZipFile(str(Path(tmpdir)) + '/proteinatlas.tsv.zip', 'r') as zipf:
        zipf.extractall(tmpdir)
    
    # get the broad expressed gene from HPA     
    hpaData = pd.read_csv(str(Path(tmpdir)) + '/proteinatlas.tsv',header=0, sep='\t')
    hpaData_broad = hpaData[hpaData['RNA tissue specificity'].isin(['Low tissue specificity','Not detected']) &\
                            hpaData['RNA tissue distribution'].isin(['Detected in many','Detected in all'])];
    hpaData_broad_gene = hpaData_broad['Ensembl'].values.tolist();
    
    broad_expressed = []
    for index, row in paralogPairs.iterrows():
        gene1 = row[0]; gene2 = row[1]; broad_exprs = 0;
        if gene1 in hpaData_broad_gene and gene2 in hpaData_broad_gene:
            broad_exprs = 1
        broad_expressed.append(broad_exprs);
    
    os.remove(str(Path(tmpdir)) + '/proteinatlas.tsv.zip')
    os.remove(str(Path(tmpdir)) + '/proteinatlas.tsv')
    return broad_expressed

if __name__ == '__main__':
    
    workdir = "/path/to/workdir/"
    
    paralog_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")

    paralog_closest = paralog_closest_check(paralogData_from_biomart=paralogData)
    paralog_chrom = paralog_chrom_check(paralogPairs=paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']])
    paralog_familysize = paralog_familysize(paralogPairs=paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']])
    

    df_merged = pd.merge(paralog_chrom, paralog_closest,  how='inner', on=['geneid','paralog_gene'])
    df_merged = pd.merge(df_merged, paralog_familysize,  how='inner', on=['geneid','paralog_gene'])
    df_merged['HPA_broad_exprs'] = paralog_HPA_exprs(paralogPairs=paralogData[['ensembl_gene_id', 'hsapiens_paralog_ensembl_gene']], 
                                                     tmpdir='/Users/cpdong/Downloads/')

    df_merged.to_csv(str(Path(workdir)) + '/Paralog_familysize_chrom_info.txt', index=False, sep='\t')
