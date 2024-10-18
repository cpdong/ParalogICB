#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 2023

@author: cpdong
"""

import argparse;
import gzip, os, shutil, time;
import pandas as pd;
import numpy as np;
from pathlib import Path
from scipy.stats import spearmanr
import urllib.request;
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser()
parser.add_argument('-c', metavar = 'cancertype', dest='cancertype', help='cancertype');
args = parser.parse_args();

'''
Use brew x86 to install x86_64 python 3.9  on M1 Chip MacOS, ref https://www.qiniu.com/qfans/qnso-70315418#comments
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
alias brew86="arch -x86_64 /usr/local/bin/brew" # add to zshrc file
brew86 install python@3.9
'''

def get_gene_pseudoKO_foldchange(data, genelist):
    #data = tcgaData.copy()
    genelist = pd.unique(paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].values.ravel('K')).tolist()
    out = pd.DataFrame()
    for gene in genelist:
        if  gene in data.index and data.loc[gene].var() != 0:
        # gene = 'ENSG00000271254'
        # geneid at index, sampleid at columns
            ctrl_sample = data.loc[gene][data.loc[gene] >= data.loc[gene].quantile(0.75)].index.tolist()
            trt_sample = data.loc[gene][data.loc[gene] <= data.loc[gene].quantile(0.25)].index.tolist()
    
            data0 = pd.DataFrame({'meanT': data[data.columns.intersection(trt_sample)].mean(axis=1),
                          'meanN': data[data.columns.intersection(ctrl_sample)].mean(axis=1)})
            out[gene] = data0['meanT'].div(data0["meanN"], axis=0)
        
        if genelist.index(gene) % 5 == 0:
            print(time.strftime("%Y-%m-%d %H:%M:%S") + ": start read!")
            print("processing: ", genelist.index(gene))
            
    out = np.log2(out)
    out.replace([np.inf, -np.inf], np.nan, inplace=True);
    return out


def calc_AUC(rank1, rank2, cutoff = 0.25): # give two rank list of genes
    cutoff = 0.25; # Pick top25 gene as true marker, here we only count the top gene not bottom ones
    
    # rank1 upon rank2
    rank1_gene_cutoff = rank1[:round(len(rank1) * cutoff)]
    rank2_gene_cutoff = rank2[:round(len(rank2) * cutoff)]
    gene_intersect = list(set(rank1_gene_cutoff) & set(rank2_gene_cutoff))

    ranking1 = [ rank2.index(x) for x in gene_intersect ]
    auc1 = (round(len(rank2) * cutoff) * len(gene_intersect)  - sum(ranking1)) / (len(rank2) * len(list(set(rank1_gene_cutoff) & set(rank2))))
    
    ranking2 = [ rank1.index(x) for x in gene_intersect ]
    auc2 = (round(len(rank1) * cutoff) * len(gene_intersect)  - sum(ranking2)) / (len(rank1) * len(list(set(rank2_gene_cutoff) & set(rank1))))
    return (auc1 + auc2)/2

def hypergeom_test(x, N, n, G, alternative='greater'):
    import scipy.stats as stats
    """
    # https://www.programcreek.com/python/example/121130/scipy.stats.hypergeom.cdf
    Parameters:
    x : int, number of `good` elements observed in the sample
    N : int, population size
    n : int, sample size
    G : int, hypothesized number of good elements in population
    alternative (default: 'greater'): {'greater', 'less', 'two-sided'}
    """
    #x=2; N=60; G=25; n=7
    pval = stats.hypergeom.sf(x-1, N, G, n)
    return pval 
    # run as hypergeom_test(x, N, n, G)

# synthetic gene pairs for survival: 0.1016/j.celrep.2019.06.067
# new scores: Mapping the landscape of synthetic lethal interactions in liver cancer

if __name__ == "__main__":
    
    workdir="/gpfs/gibbs/pi/chen_sidi/cd973/0myPrj/paralogs/data/08_feature/tcga/"
    
    paralog_file="/gpfs/gibbs/pi/chen_sidi/cd973/0myPrj/paralogs/data/07_coevolute/00_Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')
    allgenes = pd.unique(paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].values.ravel('K')).tolist()
    
    annot = pd.read_csv('https://raw.githubusercontent.com/cpdong/ParalogICB/main/data/gencode.v36_annotation.tsv',header=0,sep='\t')
    pc_genes = annot.loc[annot['biotype']=='protein_coding','Geneid'].values.tolist()

    # Go through ALL,NSCLC,RCC,BLCA,ESCA,GBM,HNSC,LIHC,SKCM,STAD
    pert_TCGA  = paralogData.copy()
    #cantypes = ['ALL', 'RCC','SKCM','BLCA','STAD','LIHC','HNSC','ESCA','GBM','NSCLC']
    #for cancer in cantypes:
    #cancer = 'GBM'
    cancer = args.cancertype;
    if cancer == 'ALL':
        tcgaFile = str(Path(workdir)) + '/TCGA-PANCAN_tumor_tpm.tsv' # note the TCGA file path
    else:
        tcgaFile = str(Path(workdir)) + '/TCGA-' + cancer + '_tumor_tpm.tsv' # note the TCGA file path

    tcgaData = pd.read_csv(tcgaFile, header=0, index_col=0, sep="\t")
    tcgaData = tcgaData.loc[pc_genes]

    fc_matrix = get_gene_pseudoKO_foldchange(tcgaData, allgenes);
    fc_matrix.to_csv(str(Path(workdir)) + '/pertKO/fc_matrix_TCGA-' + cancer + '_data.tsv', index=True, sep='\t')
       
    rhos = []; aucs = [];
    for index, row in pert_TCGA.iterrows():
        gene1 = row[0]; gene2 = row[1];
        #gene1 = 'ENSG00000000419'; gene2 = 'ENSG00000120697';
        rho = auc = None; 
        if gene1 in fc_matrix.columns and gene2 in fc_matrix.columns:
            # Spearman-test   
            df_fc = fc_matrix[[gene1,gene2]]
            df_fc.columns = ['gene1','gene2']
            df_fc.dropna(how="any", inplace=True) # drop nan values 
            (rho, pval) = spearmanr(df_fc['gene1'], df_fc['gene2'])
            
            rank1 = df_fc.sort_values('gene1', ascending=False).index.tolist()
            rank2 = df_fc.sort_values('gene2', ascending=False).index.tolist()
            # auc-recover
            auc = calc_AUC(rank1, rank2, cutoff = 0.25)
                
        if index % 5 == 0:
            print(time.strftime("%Y-%m-%d %H:%M:%S") + ": start read!")
            print("processing: ", cancer, index)
        rhos.append(rho); aucs.append(rho);
    
    pert_TCGA[cancer + '_rho'] = rhos;
    pert_TCGA[cancer + '_auc'] = aucs;

    pert_TCGA.to_csv(str(Path(workdir)) + '/pertKO/paralogGenes_perturb_TCGA-' + cancer + '_data.csv', index=False)
            
    #pert_TCGA.to_csv(str(Path(workdir)) + '/paralogGenes_perturb_TCGA_data.csv', index=False)
    #del tcgaData, pert_TCGA;
    
