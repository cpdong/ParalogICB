#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:57:44 2023

@author: cpdong
"""
import argparse, os, re;
import gzip, shutil, time;
import bs4,subprocess;
import requests
import pandas as pd;
import numpy as np;
from pathlib import Path
from functools import reduce
from scipy.stats import kendalltau
from scipy.stats import spearmanr
from scipy.stats import zscore
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.formula.api as smf
import urllib.request;

'''
Use brew x86 to install x86_64 python 3.9  on M1 Chip MacOS, ref https://www.qiniu.com/qfans/qnso-70315418#comments
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
alias brew86="arch -x86_64 /usr/local/bin/brew" # add to zshrc file
brew86 install python@3.9
'''

# Incorporate TIDE T cell dysfunction and exclude
# CD8T/NK/Macrophage propotion using GSEApy ssGSEA
# Hypoxia gene: 10.1038/s42255-019-0045-8

def TCGAbiolinks_download(project = 'ALL', dataformat = 'tpm', tissue='tumor-normal', savedir='./'):
    '''
    project like: TCGA-SKCM, MMRF-COMMPASS; if all, all 33 TCGA data set + MMRF-COMMPASS will processed
    dataformat: tpm or count
    tissue option: tumor; tumor-normal; normal;
    savedir: a directory for storage temp/final dowmloading data
    '''
    if savedir=='./':
        savedir = os.getcwd()
    os.chdir(Path(savedir))
    
    # Download TCGA data using TCGAbiolinks package in R
    os.environ["Renviron"] = '/usr/local/bin'
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    # Check the following package intalled in R    
    base = importr('base')
    TCGAbiolinks = importr('TCGAbiolinks')
    
    # Build query for ALL
    # Using Rpy2 in python: https://github.com/fcomitani/tapir wrapping mcpcounter
    if project == 'ALL':
        queryGDC = pd.DataFrame(ro.conversion.rpy2py(TCGAbiolinks.getGDCprojects()), 
                                index=list(base.colnames(TCGAbiolinks.getGDCprojects()))).T
        queryTCGA = queryGDC.loc[queryGDC['project_id'].str.contains("TCGA"), ['project_id','tumor']]
        project_ids = queryTCGA['project_id'].tolist() + ['MMRF-COMMPASS']
    elif 'TCGA' in project or 'MMRF' in project:
        project_ids = [project]
    
    assert(isinstance(project_ids, list));
    
    if os.path.exists(str(Path(savedir)) + '/MANIFEST.txt'):
        os.remove(str(Path(savedir)) + "/MANIFEST.txt")
    if os.path.exists(str(Path(savedir)) + '/GDCdata'):
        shutil.rmtree(str(Path(savedir)) + '/GDCdata')
    for project_id in project_ids:
        query = TCGAbiolinks.GDCquery(project = project_id,
                   data_category = 'Transcriptome Profiling',
                   data_type ='Gene Expression Quantification',
                   experimental_strategy='RNA-Seq',
                   workflow_type='STAR - Counts')
        TCGAbiolinks.GDCdownload(query, method = 'api', files_per_chunk = 200) # change all dot parameter to _
        queryData = ro.conversion.py2rpy(TCGAbiolinks.GDCprepare(query))
        
        getAssay = queryData.slots['assays'];
        AssayData = getAssay.slots['data'];
        listData = AssayData.slots['listData'];
        
        if dataformat == 'tpm':
            exprDat = pd.DataFrame(base.data_frame(listData[base.names(listData).index('tpm_unstrand')]),
                                 columns=base.rownames(queryData),
                                 index=base.colnames(queryData)).T
            exprDat = round(np.log2(exprDat + 1),4); # log2plus1 processing
        elif dataformat == 'count':
            exprDat = pd.DataFrame(base.data_frame(listData[base.names(listData).index('unstranded')]),
                                 columns=base.rownames(queryData),
                                 index=base.colnames(queryData)).T

        assert(isinstance(exprDat, pd.DataFrame)); # only accept tpm and count format here
        exprDat = exprDat.loc[~exprDat.index.str.contains('_PAR_')]
        exprDat.index = [x[:15] for x in exprDat.index]
        
        if 'TCGA' in project_id:
            # get only no dupicated tumor samples
            meta = pd.DataFrame(exprDat.columns,columns=['ID']);
            # remove the identical barcode as TCGA-MJ-A850-01A TCGA-MJ-A850-01B or TCGA-MJ-A850-01C
            meta['sID'] = [x[:15] for x in meta.ID]; meta['sID_rank'] = [ord(x[15:16].upper())-64 for x in meta.ID]
            meta = meta.loc[meta.groupby('sID').sID_rank.idxmin()]
            meta['tissue'] = [x[13:15] for x in meta.ID];
            #meta.tissue.value_counts()
            normal_samples = meta.loc[meta['tissue'] == '11','ID'].tolist(); # keep for future option
            meta = meta.loc[meta['tissue'] != '11'] # remove normal samples
            # remove the identical barcode as TCGA-MJ-A850-01 TCGA-MJ-A850-06, get in a 01 02 03 06 order
            meta['pID'] = [x[:12] for x in meta.ID];
            meta['pID_rank'] = [int(x[13:15]) for x in meta.ID]
            meta = meta.loc[meta.groupby('pID').pID_rank.idxmin()]
            tumor_samples = meta['ID'].tolist();
            
            if tissue == 'tumor':
                exprDat2 = exprDat[exprDat.columns.intersection(tumor_samples)]
            elif tissue == 'normal':
                exprDat2 = exprDat[exprDat.columns.intersection(normal_samples)]
            elif tissue == 'tumor-normal':
                exprDat2 = exprDat[exprDat.columns.intersection(tumor_samples + normal_samples)]
            else:
                print('TCGAbiolinks: Check your tissue setting!')
                
        elif project_id == 'MMRF-COMMPASS':
            meta = pd.DataFrame(exprDat.columns,columns=['ID']);
            meta['tissue'] = [x[10:14] for x in meta.ID];
            normal_samples = []; # no normal sample in MMRF
            tumor_samples = meta.loc[meta['tissue'] == '1_BM','ID'].tolist(); # tumor refer to primary NDMM
            # secondary = meta.loc[meta['tissue'] == '2_BM','ID'].tolist();

            if tissue == 'tumor':
                exprDat2 = exprDat[exprDat.columns.intersection(tumor_samples)]
            elif tissue == 'tumor-normal': # here output all 859 samples expression profile
                exprDat2 = exprDat
            else:
                print('TCGAbiolinks: Check your tissue setting for MMRF! No normal sample in MMRF.')
        
        exprDat2.index.rename('Geneid', inplace=True)
        exprDat2.to_csv(str(Path(savedir)) + '/' + project_id + '_' + tissue + '_' + dataformat + '.tsv', index=True, sep='\t');
        
        os.remove(str(Path(savedir)) + "/MANIFEST.txt")
        shutil.rmtree(str(Path(savedir)) + "/GDCdata")

#TCGAbiolinks_download(project = 'TCGA-CHOL', dataformat = 'tpm',  tissue='tumor', savedir='/Users/cpdong/Downloads/tt/')     

def gsvar(data, gsets): # wrapper orignial gsva function from R
    '''
    R and python at the same conda environment
    data: pandas df
    gsets in a dictionary: gene_sets={'A':['gene1', 'gene2',...]}
    '''
    # Running GSVA/ssGSEA using GSVA package in R
    os.environ["Renviron"] = '/usr/local/bin'
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    # Check the following package intalled in R
    base = importr('base')
    GSVA = importr('GSVA')
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        tabr = ro.conversion.py2rpy(data)

    gsets_list = ro.ListVector([(str(k), x) for k, x in gsets.items()])
    
    ssGSEA_Score = GSVA.gsva(expr = base.data_matrix(tabr),
                       gset_idx_list = gsets_list,
                       method = 'ssgsea',
                       kcdf = 'Gaussian') # default abs.ranking=FALSE

    ssgsea_df = pd.DataFrame(base.data_frame(ssGSEA_Score),
                         columns=base.rownames(ssGSEA_Score),
                         index=base.colnames(ssGSEA_Score))
    return ssgsea_df

def hypoxia_score(data):
    '''
    The 15 hypoxia gene from 10.1038/s42255-019-0045-8
    '''
    hypoxia_df = pd.read_csv("https://raw.githubusercontent.com/cpdong/public/master/Hypoxia_GeneSignature_pmid31984309.csv", header=0)
    gsets={'Hypoxia': hypoxia_df['Geneid'].values.tolist()}
    
    hypoxia_df = gsvar(data, gsets)
    
    return hypoxia_df

def immune_score(data):
    '''
    The 17 immune cell signature from pmid28052254
    Here we only keep: CD8T/NK/Macrophage::CD8 T cell activated; NK activated; 
    '''
    immune_gene = pd.read_csv("https://raw.githubusercontent.com/cpdong/public/master/ImmCellSig_TNKMacro_signatures.csv", header=0)
    immune_dict = {k: list(v) for k, v in immune_gene.groupby('name')['ensembl_gene_id']}
    # immune_dict = immune_gene.groupby('name')['ensembl_gene_id'].apply(list).to_dict()

    immune_score = gsvar(data, immune_dict)
    
    return immune_score

def tide_score(data, cancertype='Other'):
    '''
    https://github.com/jingxinfu/TIDEpy
    Here we only keep: T cell Dysfunction/Exclusion
    '''
    from tidepy.pred import TIDE
    result = TIDE(data, cancer=cancertype, pretreat=False, vthres=0.)
    
    return result[['Dysfunction','Exclusion']]

def quantile_normalize(df):
    # Ref https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/
    # Step 1: Order values in each column
    df_sorted = pd.DataFrame(np.sort(df.values, axis=0), 
                         index=df.index, 
                         columns=df.columns)
    # Step2: Compute Row Means
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    # Step3: Use Average Values to each sample in the original order
    df_quantiled =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_quantiled)

def paralog_TME_tau(genePair_df, geneExpData):
    #genePair_df = paralogData.copy()
    #geneExpData = mixture.copy()
    genePair_df.columns = ['gene1','gene2']
    # Get hypoxia/TNKMacro/tide
    hypoxia_data = hypoxia_score(geneExpData.copy())
    immune_data = immune_score(geneExpData.copy())
    tide_data =  tide_score(geneExpData.copy())
    score0 = pd.merge(hypoxia_data, immune_data, left_index=True, right_index=True, how='inner')
    TME_scores = pd.merge(score0, tide_data, left_index=True, right_index=True, how='inner')

    exprT = geneExpData.T; del geneExpData # save memory
    
    for item in TME_scores.columns.tolist():
        tau_list = [];
        for index, row in genePair_df.iterrows():
            tau = None;
            gene1 = row[0]; gene2 = row[1];
            if gene1 in exprT.columns and gene2 in exprT.columns:
                dat1 = exprT[[gene1,gene2]]
                dat1 = (dat1 >= dat1.median()).astype('int')
                dat1['group'] = dat1[gene1] + dat1[gene2] + 1
                dat1.index.rename('Geneid', inplace=True)
                dat2 = pd.merge(dat1, TME_scores,left_index=True,right_index=True, how='inner')
                # kendall-test4
                tau, pval = kendalltau(dat2['group'], dat2[item])
            tau_list.append(tau)
    
            if index % 1000 == 0:
                print("processing: ", index)
        genePair_df[item] = tau_list;
    return genePair_df


if __name__ == "__main__":
    
    workdir="/Users/cpdong/Dropbox/project/paralogs/data/08_features"
    paralog_file="/Users/cpdong/Dropbox/project/paralogs/data/07_coevolute/00_Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

    # Download the TCGA and MMRF datasets
    TCGAbiolinks_download(project = 'ALL', dataformat='tpm', tissue='tumor', savedir='/Users/cpdong/Downloads/tt/')
    
    '''
    # Create merged dataset of TCGA-MMRF-PANCAN /TCGA-NSCLC/TCGA-RCC
    filedir='/Users/cpdong/Downloads/tt/'
    files = [x for x in os.listdir(filedir) if '_tumor_tpm.tsv' in x]
    d0 = pd.read_csv(filedir + '/' + files[0], header=0, sep='\t')
    for f in files[1:]:
        print(f)
        d1 = pd.read_csv(filedir + '/' + f, header=0, sep='\t')
        d0 = pd.merge(d0, d1, on='Geneid', how='outer')
    d0 = d0.set_index('Geneid')
    df_tcga = quantile_normalize(d0)
    df_tcga.index.rename('Geneid', inplace=True)
    df_tcga.to_csv(filedir + '/TCGA-PANCAN_tumor_tpm.tsv', index=True, sep='\t');
    del df_tcga
    
    # Create merged dataset of TCGA-NSCLC
    filedir='/Users/cpdong/Downloads/tt/'
    files = [x for x in os.listdir(filedir) if 'TCGA-LU' in x]
    d0 = pd.read_csv(filedir + '/' + files[0], header=0, sep='\t')
    for f in files[1:]:
        print(f)
        d1 = pd.read_csv(filedir + '/' + f, header=0, sep='\t')
        d0 = pd.merge(d0, d1, on='Geneid', how='outer')
    d0 = d0.set_index('Geneid')
    tcga_nsclc = quantile_normalize(d0)
    tcga_nsclc.index.rename('Geneid', inplace=True)
    tcga_nsclc.to_csv(filedir + '/TCGA-NSCLC_tumor_tpm.tsv', index=True, sep='\t');
    del tcga_nsclc
    
    # Create merged dataset of TCGA-RCC
    filedir='/Users/cpdong/Downloads/tt/'
    files = [x for x in os.listdir(filedir) if 'TCGA-KI' in x]
    d0 = pd.read_csv(filedir + '/' + files[0], header=0, sep='\t')
    for f in files[1:]:
        print(f)
        d1 = pd.read_csv(filedir + '/' + f, header=0, sep='\t')
        d0 = pd.merge(d0, d1, on='Geneid', how='outer')
    d0 = d0.set_index('Geneid')
    tcga_rcc = quantile_normalize(d0)
    tcga_rcc.index.rename('Geneid', inplace=True)
    tcga_rcc.to_csv(filedir + '/TCGA-RCC_tumor_tpm.tsv', index=True, sep='\t');
    del tcga_rcc
    '''
    
    # Go through TCGA- ALL,ALLMMRF,NSCLC,RCC,BLCA,ESCA,GBM,HNSC,LIHC,SKCM,STAD
    cantypes = ['RCC','SKCM','BLCA','STAD','LIHC','HNSC','ESCA'] #PANCAN/GBM/NSCLC
    for cancer in cantypes:
        #cancer = 'RCC'
        file = '/Users/cpdong/Downloads/tcga/TCGA-' + cancer + '_tumor_tpm.tsv'
        mixture = pd.read_csv(file, header=0, index_col=0, sep='\t')
        # mixture.index = [x[:15] for x in mixture.index]

        paras = paralog_TME_tau(paralogData.copy(), mixture);
        paras.to_csv('/Users/cpdong/Downloads/tcga/TCGA-' + cancer + '_TME_kendalltest.tsv', index=False, sep='\t');
        
