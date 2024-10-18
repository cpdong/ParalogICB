#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:57:44 2023

@author: cpdong
"""

import argparse, glob, os, re;
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

# synthetic gene pairs for survival: 0.1016/j.celrep.2019.06.067
# new scores: Mapping the landscape of synthetic lethal interactions in liver cancer


def TCGAclinical_download(project = 'ALL', savedir='./'):
    '''
    project like: TCGA-SKCM, MMRF-COMMPASS; if all, all 33 TCGA data set + MMRF-COMMPASS will processed
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
    for project_id in project_ids:
        #project_id = 'MMRF-COMMPASS'
        print(project_id)
        #savedir='/Users/cpdong/Downloads/tcga/tttt/'
        os_url = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/' + project_id + '.survival.tsv'
        urllib.request.urlretrieve(os_url, str(Path(savedir)) + '/' + project_id + '.surv.tsv'); # download surv data
        surv = pd.read_csv(str(Path(savedir)) + '/' + project_id + '.surv.tsv', header=0, sep='\t', quotechar='"')[['sample','OS.time','OS']]
        
        if project_id == 'MMRF-COMMPASS':
            phen_url = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/' + project_id + '.Xena_phenotype.tsv'
            urllib.request.urlretrieve(phen_url, str(Path(savedir)) + '/' + project_id + '.phenotype.tsv'); # download surv data
            pheno = pd.read_csv(str(Path(savedir)) + '/' + project_id + '.phenotype.tsv', header=0, sep='\t', quotechar='"')
            pheno = pheno[['samples.submitter_id','diagnoses.age_at_diagnosis','demographic.gender','demographic.race']]
            pheno.columns = ['sample','age','sex','race']
            os.remove(str(Path(savedir)) + '/' + project_id + '.phenotype.tsv')
        else:
            phen_url = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/' + project_id + '.GDC_phenotype.tsv.gz'
            urllib.request.urlretrieve(phen_url, str(Path(savedir)) + '/' + project_id + '.phenotype.tsv.gz'); # download surv data
            pheno = pd.read_csv(str(Path(savedir)) + '/' + project_id + '.phenotype.tsv.gz', compression='gzip', header=0, sep='\t', quotechar='"')
            pheno = pheno[['submitter_id.samples','age_at_initial_pathologic_diagnosis','gender.demographic','race.demographic']]
            pheno.columns = ['sample','age','sex','race']
            os.remove(str(Path(savedir)) + '/' + project_id + '.phenotype.tsv.gz')
        
        surv_pheno = pd.merge(surv, pheno, on='sample', how='outer')
        surv_pheno.to_csv(str(Path(savedir)) + '/' + project_id + '.survdata.tsv', index=False, sep='\t');
        if project_id == 'MMRF-COMMPASS':
            surv_pheno['age'] = round(surv_pheno['age'] / 365).astype('Int64')

        surv_pheno.to_csv(str(Path(savedir)) + '/' + project_id + '.survdata.tsv', index=False, sep='\t');
        os.remove(str(Path(savedir)) + '/' + project_id + '.surv.tsv')
        
    # Gather all phenotype
    all_files = glob.glob(str(Path(savedir)) + '/*.survdata.tsv')
    filelist = []
    for filename in all_files:
        dat = pd.read_csv(filename, index_col=None, header=0, sep='\t')
        filelist.append(dat)
    clinicALL = pd.concat(filelist, axis=0, ignore_index=True)
    clinicALL.to_csv(str(Path(savedir)) + '/TCGA-PANCAN.survdata.tsv', index=False, sep='\t');
    
    # Gather LUAD LUSC 
    lung_files = [x for x in all_files if 'LUAD' in x or 'LUSC' in x]
    filelist = []
    for filename in lung_files:
        dat = pd.read_csv(filename, index_col=None, header=0, sep='\t')
        filelist.append(dat)
    clinicLUNG = pd.concat(filelist, axis=0, ignore_index=True)
    clinicLUNG.to_csv(str(Path(savedir)) + '/TCGA-NSCLC.survdata.tsv', index=False, sep='\t');
    
    # Gather KIRP KIRC KICH 
    rcc_files = [x for x in all_files if 'KIRC' in x or 'KIRP' in x or 'KICH' in x]
    filelist = []
    for filename in rcc_files:
        dat = pd.read_csv(filename, index_col=None, header=0, sep='\t')
        filelist.append(dat)
    clinicRCC = pd.concat(filelist, axis=0, ignore_index=True)
    clinicRCC.to_csv(str(Path(savedir)) + '/TCGA-RCC.survdata.tsv', index=False, sep='\t');
        
#TCGAclinical_download(project = 'ALL', savedir='/Users/cpdong/Downloads/tcga/clinic/')

def get_paired_surv(geneData, phenoData, genedPairs):
    #genedPairs = paralogData.copy()
    #geneData = pd.read_csv('/Users/cpdong/Downloads/tcga/TCGA-BLCA_tumor_tpm.tsv', index_col=0, header=0, sep='\t')
    geneData.columns = [x[:16] for x in geneData.columns]
    #phenoData = pd.read_csv('/Users/cpdong/Downloads/tcga/clinic/TCGA-PANCAN.survdata.tsv', index_col=0, header=0, sep='\t')
    phenoData['race'] = phenoData['race'].map({'white': 'white',
                                               'asian': 'asian',
                                               'black or african american': 'black',
                                               'not reported': None,
                                               'american indian or alaska native': 'other',
                                               'native hawaiian or other pacific islander': 'other',
                                               'not allowed to collect': None,
                                               'other': 'other'})
    phenoData['sex'] = phenoData['sex'].map({'female': 'female', 'male': 'male', 'not reported': None})
    #phenoData.to_csv(str(Path(savedir)) + '/TCGA-ALL.survdata.tsv', index=False, sep='\t');
    
    data = pd.merge(geneData.T, phenoData, left_index=True, right_index=True, how='inner')
    
    # Running survival using survival package in R
    os.environ["Renviron"] = '/usr/local/bin'
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.conversion import localconverter
    import warnings
    from pandas.errors import SettingWithCopyWarning
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
    
    base = importr('base')
    survival = importr('survival')
    survscore = []
    for index, row in genedPairs.iterrows():
        gene1 = row[0]; gene2 = row[1]; score = None;
        #gene1='ENSG00000000419'; gene2='ENSG00000000457'
        if gene1 in data.columns and gene2 in data.columns:
            sdata = data[[gene1, gene2, 'OS.time', 'OS', 'age', 'sex', 'race']]
            sdata['group'] = np.where((sdata[gene1] >= sdata[gene1].median()) & (sdata[gene2] >= sdata[gene2].median()), 3,
                                      np.where((sdata[gene1] < sdata[gene1].median()) & (sdata[gene2] < sdata[gene2].median()), 1, 2))
            if len(sdata['group'].unique()) == 3: # get correct grouping, skip all 0 genes.
                sdata_d1 = sdata[sdata['group'] != 1];
                sdata_d2 = sdata[sdata['group'] != 3];
                with localconverter(ro.default_converter + pandas2ri.converter):
                    survdata = ro.conversion.py2rpy(sdata)
                    surv_d1 = ro.conversion.py2rpy(sdata_d1)
                    surv_d2 = ro.conversion.py2rpy(sdata_d2)
                '''
                # Build command
                robjects.globalenv['r_df'] = survdata
                cmd = ('coxfit = coxph(Surv(r_df${time}, r_df${event}) ~ {incols}, '
                       + 'data=r_df, model=TRUE)').format(time='OS.time',
                                                          event='OS',
                                                          incols='group + age + sex')
                coxfit = r(cmd)
                print(r('summary(coxfit)$coefficients'))
                '''
                # chceked with R coxph function
                coxph = pd.DataFrame(ro.conversion.py2rpy(base.summary(survival.coxph(Formula("Surv(OS.time, OS) ~ group + age + sex"), data = survdata))[6]))
                coxph3v2= pd.DataFrame(ro.conversion.py2rpy(base.summary(survival.coxph(Formula("Surv(OS.time, OS) ~ group"), data = surv_d1))[6]))
                coxph2v1= pd.DataFrame(ro.conversion.py2rpy(base.summary(survival.coxph(Formula("Surv(OS.time, OS) ~ group"), data = surv_d2))[6]))
        
                cox_p = coxph[0][(len(coxph[0])/5)*4]
                sign = np.sign((coxph3v2[0][1] -1) * (coxph2v1[0][1] -1))
                score = round(-sign * np.log10(cox_p),4)
           
        survscore.append(score)
        
        if index % 1000 == 0:
            print('processing:', index)
        
    return survscore
    
#TCGAclinical_download(project = 'ALL', savedir='/Users/cpdong/Downloads/tcga/tttt/')

if __name__ == "__main__":
    
    workdir="/path/to/workdir/"
    paralog_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

    TCGAclinical_download(project = 'ALL', savedir='/path/to/workdir/clinic/')
    # Go through TCGA- ALL,ALLMMRF,NSCLC,RCC,BLCA,ESCA,GBM,HNSC,LIHC,SKCM,STAD
    
    paralogDF = paralogData.copy()
    cantypes = ['PANCAN','NSCLC','RCC','GBM','SKCM','BLCA','STAD','LIHC','HNSC','ESCA']
    for cancer in cantypes:
        #cancer = 'BLCA'
        genedPairs = paralogData.copy()
        geneData = pd.read_csv(str(Path(workdir)) + '/TCGA-' + cancer + '_tumor_tpm.tsv', index_col=0, header=0, sep='\t')
        phenoData = pd.read_csv(str(Path(workdir)) + '/clinic/TCGA-' + cancer  + '.survdata.tsv', index_col=0, header=0, sep='\t')
        
        survscores = get_paired_surv(geneData, phenoData, genedPairs)
        paralogDF['cancer' + '_surv'] = survscores
        paralogDF.to_csv(str(Path(workdir)) + '/paralogs_survscores.tsv', index=False, sep='\t');
