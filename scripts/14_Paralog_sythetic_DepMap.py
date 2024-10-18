#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:57:44 2023

@author: cpdong
"""
import argparse, os, re;
import gzip, shutil, time;
import bs4
import requests
import pandas as pd;
import numpy as np;
from pathlib import Path
from scipy.stats import spearmanr
from scipy.stats import zscore
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.formula.api as smf
import urllib.request;

workdir="/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/08_features"
cancerType="SKCM"

paralog_file="/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/07_coevolute/00_Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

DepMap_gene_IDmapping="https://raw.githubusercontent.com/cpdong/public/master/DepMap_IDmapping.tsv"
DepMap_gene_IDmap = pd.read_csv(DepMap_gene_IDmapping,header=0,sep="\t")
DepMap_gene_IDmap = DepMap_gene_IDmap[DepMap_gene_IDmap['Geneid'].notnull()]

'''
# DepMap Public 23Q2 Files Download automatically
DEPMAP_API_URL = "https://depmap.org/portal/download/api/downloads"

filelist = ['CRISPRGeneDependency.csv',
            'CRISPRGeneEffect.csv',
            'ScreenGeneDependency.csv',
            'ScreenGeneEffect.csv',
            'OmicsExpressionProteinCodingGenesTPMLogp1.csv',
            'OmicsCNGene.csv',
            'OmicsSomaticMutations.csv',
            'Model.csv']

depmap_datadir = str(Path(workdir)) + "/DepMap_data"
if not os.path.exists(depmap_datadir):
    os.makedirs(depmap_datadir)
    
# Ref https://github.com/cthoyt/depmap_downloader
depmap_json_table = requests.get(DEPMAP_API_URL).json()
latest = next( release for release in depmap_json_table["releaseData"] if release["isLatest"] )
#version = latest["releaseName"] # Format as 'DepMap Public 22Q4', you can change here
version = 'DepMap Public 23Q2'


for fname in filelist:
    for download in depmap_json_table["table"]:
        if download["fileName"] == fname and download["releaseName"] == version:
            real_url = download["downloadUrl"]
            full_name = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_" + fname
            urllib.request.urlretrieve(real_url, full_name)
'''
depmap_datadir = str(Path(workdir)) + "/DepMap_data"
version = 'DepMap Public 23Q2'

CRISPRGeneDependency = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_CRISPRGeneDependency.csv"
CRISPRGeneEffect = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_CRISPRGeneEffect.csv"
ScreenGeneDependency = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_ScreenGeneDependency.csv"
ScreenGeneEffect = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_ScreenGeneEffect.csv"
OmicsCNV = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsCNGene.csv"
OmicsExpr = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsExpressionProteinCodingGenesTPMLogp1.csv"
OmicsSNV = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsSomaticMutations.csv"
cell_annotation = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_Model.csv"

# DepMap cell match TCGA abbrev code
cellAnnot = pd.read_csv(cell_annotation, header=0, index_col=0)
DepMap_cell_TCGAmapping="https://raw.githubusercontent.com/cpdong/public/master/DepMap_TCGA_cell_annotation.tsv"
DepMap_cell_TCGAmap = pd.read_csv(DepMap_cell_TCGAmapping, header=0, sep="\t")
DepMap_cell_TCGAmap = DepMap_cell_TCGAmap[~DepMap_cell_TCGAmap['TCGAcode'].isin(['COAD','READ','LUAD','LUSC','KIRC','KIRP'])] # anno as COADREAD,NSCLC,RCC
cell_reAnnot = pd.DataFrame(columns=list(cellAnnot.columns) + ['TCGAcode'])
for index, row in DepMap_cell_TCGAmap.iterrows():
    cantype = DepMap_cell_TCGAmap['TCGAcode'][index]
    mydict = DepMap_cell_TCGAmap.loc[index,].to_dict()
    mydict = {k: v for k, v in mydict.items() if k != 'TCGAcode' and not pd.isna(v) }

    cellAnnot0 = cellAnnot
    for filters, values in mydict.items():
        cellAnnot0 = cellAnnot0[cellAnnot0[filters].isin(values.split('|'))]
    cellAnnot0['TCGAcode'] = cantype
    cell_reAnnot = pd.concat([cell_reAnnot, cellAnnot0])
cell_reAnnot.to_csv(str(Path(depmap_datadir)) + '/DepMap_TCGAcell_annotation.txt', index=True, sep='\t')


# Using Archilies CRISPR data
# Essential gene = CERES score < -0.6, as in (Meyers et al., 2017)
# Synthetic gene method as 10.1016/j.cels.2021.08.006

# CRISPR or Screen

# CRISPRGeneEffect parsing
CRISPRGeneEffs = pd.read_csv(CRISPRGeneEffect,header=0,index_col=0)
CRISPRGeneEffs = CRISPRGeneEffs[CRISPRGeneEffs.columns.intersection(DepMap_gene_IDmap['ID'])].dropna(axis=1, how='any') # drop column with NaN, esp for ScreenGeneEff data
CRISPRGeneEffs = CRISPRGeneEffs.rename(columns=DepMap_gene_IDmap.set_index('ID')['Geneid']) # replace ID to ensemblid
CRISPRGeneEffs = CRISPRGeneEffs[CRISPRGeneEffs.columns.intersection(list(pd.unique(paralogData.values.ravel('K'))))] # get all paralog genes only
# only 27 cantype match with TCGA code
CRISPRGeneEffs_tcga = CRISPRGeneEffs[CRISPRGeneEffs.index.isin(cell_reAnnot.index)]


# Func-1: essential genes
'''
# CRISPRGeneEff.to_csv(str(Path(depmap_datadir)) + '/DepMap_SKCM_eff.txt', index=True, sep='\t')

# Essential gene function from shinyDepMap:10.7554/eLife.57116; 10.1038/s41598-022-16940-7
# The above method can be used but when used with GeneEffect value, seems incorrect
# OR with a specific cutoff criteria: Essential in ≥ 10% cell lines(10.1101/2022.06.12.495566)
# Gene essentiality was determined by the mean Chronos score across 34 MM cell lines(10.1101/2023.04.04.535554)
# Gene to be considered essential at least 1% of all cell lines(CHRONOS score < −0.6, doi.org/10.1101/2023.03.02.530664)
# Here we try to use chronos -0.6, 10% cell lines as well
'''
def paralog_essential(paralogData,ScreenGeneEffData, ess_cutoff, perc_cutoff):
    # paralogData - a two cols DF, col1 paralog1, col2 paralog2
    # ScreenGeneEffData download from DepMap, not dependency data
    # ess_cutoff, CERES and Chronos cutoff -0.6
    # perc_cutoff, define essital for percentage cutoff of cell with score less than ess_cutoff
    # ess_cutoff = -0.6; perc_cutoff=0.1;
    # single gene functions
    # es_df = CRISPRGeneEff[CRISPRGeneEff < ess_cutoff].count(axis=0)/CRISPRGeneEff.shape[0]
    # es_gene = es_df[es_df>0.1].index
    
    # paired gene functions, either one <-0.6 define essential
    essentials = [];
    for index, row in paralogData.iterrows():
        essential_eval = None
        #gene1 = row['ensembl_gene_id'];
        #gene2 = row['hsapiens_paralog_ensembl_gene'];
        gene1 = row[0]; gene2 = row[1];
        if any(x in [gene1, gene2] for x in list(ScreenGeneEffData.columns)):
            parag_check = ScreenGeneEffData[ScreenGeneEffData.columns.intersection([gene1,gene2])]
            parag_check_min = parag_check.min(axis=1)
            ess_cnts = len(parag_check_min[parag_check_min < ess_cutoff])
            if ess_cnts/ScreenGeneEffData.shape[0] > perc_cutoff:
                essential_eval = 1
            else:
                essential_eval = 0
        essentials.append(essential_eval)
        
    return essentials

Result_Essential = paralogData.copy()
Result_Essential['cell_all_essential'] = paralog_essential(paralogData, CRISPRGeneEffs, ess_cutoff=-0.6, perc_cutoff=0.01)
Result_Essential['cell_tcga_essential'] = paralog_essential(paralogData, CRISPRGeneEffs_tcga, ess_cutoff=-0.6, perc_cutoff=0.01)
# Go to find essentials for per Cancer Type
for cancerType in ['BLCA','ESCA','GBM','HNSC','RCC','LIHC','NSCLC','SKCM','STAD','BRCA','COADREAD']:
    print("processing essential genes for", cancerType)
    cell_canType = cell_reAnnot[cell_reAnnot['TCGAcode'] == cancerType]
    CRISPRGeneEff = CRISPRGeneEffs[CRISPRGeneEffs.index.isin(cell_canType.index)]
    Result_Essential['cell_' + cancerType + '_essential'] = paralog_essential(paralogData,CRISPRGeneEff,ess_cutoff=-0.6, perc_cutoff=0.1)

Result_Essential.to_csv(str(Path(depmap_datadir)) + '/DepMap_paralogGenes_essential_data.csv', index=False)

# Func-2: sythetic genes analysis
'''
Step-1: prepare the loss gene matrix 
'''

OmicsCNV = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsCNGene.csv"
OmicsExpr = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsExpressionProteinCodingGenesTPMLogp1.csv"
OmicsSNV = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsSomaticMutations.csv"

# collect gene loss from gene expression data
exprDat = pd.read_csv(OmicsExpr,header=0,index_col=0)
exprDat = exprDat[exprDat.columns.intersection(DepMap_gene_IDmap['ID'])]
exprDat = exprDat.rename(columns=DepMap_gene_IDmap.set_index('ID')['Geneid']) # replace ID to ensemblid
lossExpTpm = pd.DataFrame(data=np.where(exprDat < 1, 1, 0), index=exprDat.index, columns=exprDat.columns) # log2TPMp1<1 as loss
# zscore - step
# exprDat = exprDat.T.apply(zscore).T # zscore process for each sample
exprDat = exprDat.apply(zscore) # zscore process for each gene
lossExpZscore = pd.DataFrame(data=np.where(exprDat < -4, 1, 0), index=exprDat.index, columns=exprDat.columns) # zscore<-4 as loss
# combine geneTPM and z-scores
lossExpSum = lossExpTpm.add(lossExpZscore, fill_value=0) # combine gene expression logTPM and zscore

# collect gene loss from copy number data
cnvDat0 = pd.read_csv(OmicsCNV,header=0,index_col=0)
cnvDat = np.log2(np.exp2(cnvDat0) - 1) # convern log2(x+1) to log2(x)
cnvDat = cnvDat[cnvDat.columns.intersection(DepMap_gene_IDmap['ID'])]
cnvDat = cnvDat.rename(columns=DepMap_gene_IDmap.set_index('ID')['Geneid']) # replace ID to ensemblid
# cnvDat = cnvDat[cnvDat.index.isin(cell_canType['ModelID'])]
lossCNV = pd.DataFrame(data=np.where(cnvDat < -1.28, 1, 0), index=cnvDat.index, columns=cnvDat.columns) # log2CNV<-1.28 as loss

comm_gene = list(set(cnvDat.columns) & set(lossExpSum.columns))
comm_cell = list(set(cnvDat.index) & set(lossExpSum.index))
lossCNV2 = lossCNV.loc[comm_cell][comm_gene]
lossExpSum2 = lossExpSum.loc[comm_cell][comm_gene]
lossExpCnv = lossExpSum2.add(lossCNV2, fill_value=0) # combine gene expression and CNV loss informaton

# collect gene loss from mutation info
snvDat = pd.read_csv(OmicsSNV,header=0)
snvDat.rename(columns={'ModelID':'DepMap_ID'}, inplace=True); # some version colname changed
snvDat = snvDat[snvDat['VariantInfo'].isin(['SPLICE_SITE','NONSENSE','FRAME_SHIFT_DEL','FRAME_SHIFT_INS'])] # get only LoF mutation
snvDat = snvDat[snvDat.columns.intersection(['DepMap_ID','EntrezGeneID'])].drop_duplicates(keep='first')
snvDat = snvDat[(snvDat['EntrezGeneID'].notna()) & (snvDat['EntrezGeneID'] !='Unknown')]
snvDat.EntrezGeneID = snvDat.EntrezGeneID.astype(float).astype(int)
# convert EntrezGeneID to ensembl id
anno2 = DepMap_gene_IDmap[['Geneid','entrezgene_id']].groupby(['Geneid'], as_index=False)['entrezgene_id'].min()
snvDat2 = pd.merge(snvDat, anno2, left_on='EntrezGeneID', right_on='entrezgene_id', how='inner')
lossSnv = pd.DataFrame(data=np.where(lossExpCnv > 0, 0, 0), index=lossExpCnv.index, columns=lossExpCnv.columns)
for index, row in snvDat2.iterrows():
    #cellid = row['DepMap_ID']; #geneid = row['Geneid'];
    if row['DepMap_ID'] in lossSnv.index and row['Geneid'] in lossSnv.columns:
        lossSnv.loc[row['DepMap_ID'], row['Geneid']] = 1

lossExpCnvSnv = lossExpCnv.add(lossSnv, fill_value=0) # combine gene expression, CNV and mutation loss informaton
lossExpCnvSnv2 = pd.DataFrame(data=np.where(lossExpCnvSnv > 0, 1, 0), index=lossExpCnvSnv.index, columns=lossExpCnvSnv.columns).astype('Int64')
lossExpCnvSnv2 = lossExpCnvSnv2[lossExpCnvSnv2.columns.intersection(list(pd.unique(paralogData.values.ravel('K'))))] # get all paralog genes only
lossExpCnvSnv2 = lossExpCnvSnv2[lossExpCnvSnv2.index.isin(CRISPRGeneEffs.index)] # find overlap cell with CRIPSR and gene loss information for fitting models

lossExpCnvSnv2.to_csv(str(Path(depmap_datadir)) + '/DepMap_geneLoss_data.csv', index=True)

# Check all loss statistics for each round
lossExpTpm.stack().value_counts() # check only tpm loss cases
lossExpZscore.stack().value_counts() # check only gene zscore loss cases
lossCNV.stack().value_counts() # check only CNV loss cases
lossSnv.stack().value_counts() # check only mutatiopn loss cases
lossExpCnv.stack().value_counts()  # check Exp + CNV loss cases
lossExpCnvSnv2.stack().value_counts()  # check final loss cases


# Single gene essential for sythetic calling
def single_gene_essential(genelist, ScreenGeneEffData, ess_cutoff, perc_cutoff):
    essentials = []; chronos_avg = []
    for gene in genelist:
        essential_eval = None; chronos = None;
        if gene in list(ScreenGeneEffData.columns):
            chronos = ScreenGeneEffData[gene].replace('NaN', np.nan).mean(skipna=True)
            geneEff = ScreenGeneEffData[gene]
            ess_cnts = len(geneEff[geneEff < ess_cutoff])
            if ess_cnts/ScreenGeneEffData.shape[0] > perc_cutoff:
                essential_eval = 1
            else:
                essential_eval = 0
        chronos_avg.append(chronos)
        essentials.append(essential_eval)
    return essentials, chronos_avg

# calc essential for each gene
genelist = list(pd.unique(paralogData.values.ravel('K'))) # get all gene in the paralogData
geneEssential = pd.DataFrame(genelist, columns=['Geneid'])
geneEssential['cell_all_essential'], geneEssential['cell_all_avgChronos'] = single_gene_essential(genelist, CRISPRGeneEffs, ess_cutoff=-0.6, perc_cutoff=0.01)
geneEssential['cell_tcga_essential'], geneEssential['cell_tcga_avgChronos'] = single_gene_essential(genelist, CRISPRGeneEffs_tcga, ess_cutoff=-0.6, perc_cutoff=0.01)
# Go to find essentials for per Cancer Type
for cancerType in ['BLCA','ESCA','GBM','HNSC','RCC','LIHC','NSCLC','SKCM','STAD','BRCA','COADREAD']:
    print("processing single gene essential genes for", cancerType)
    cell_canType = cell_reAnnot[cell_reAnnot['TCGAcode'] == cancerType]
    CRISPRGeneEff = CRISPRGeneEffs[CRISPRGeneEffs.index.isin(cell_canType.index)]
    geneEssential['cell_' + cancerType + '_essential'], geneEssential['cell_' + cancerType + '_avgChronos'] = single_gene_essential(genelist,CRISPRGeneEff,ess_cutoff=-0.6, perc_cutoff=0.1)

geneEssential.to_csv(str(Path(depmap_datadir)) + '/DepMap_singleGene_essential_avgChronos_data.csv', index=False)

def paralog_synthetic(paralogPair, dependData, genelossData, geneEss_dict, cellAnnot):
    #paralogPair = paralogData
    #dependData = CRISPRGeneDependency
    #genelossData = lossExpCnvSnv2
    #geneEss_dict = dict(zip(geneEssential['Geneid'], geneEssential['cell_all_essential']))
    #cellAnnot = cellAnnot
    #cellAnnot['celltype'] = cellAnnot['OncotreeLineage']
    #celltypeDat = cellAnnot[['OncotreeLineage']]; celltypeDat.columns = ['celltype']
    #cellAnnot = celltypeDat
    
    # SL criteria: glm coef < 0 and adjust p<0.1 as Cell syst paper
    #paralogPair = paralogData
    para1 = paralogPair.copy(); para1.columns = ['gene1','gene2']
    para2 = para1[['gene2','gene1']]; para2.columns = ['gene1','gene2']
    paralogPairs = pd.concat([para1, para2]).drop_duplicates(keep='first') # switch A1 A2, rbind and remove duplicated rows
    paralogPairs.reset_index(inplace=True, drop=True) # reindexing
    
    paralogPairs['coef'] = None;paralogPairs['pval'] = None;paralogPairs['synthetic_p005'] = None;
    paralogPairs['FDR'] = None;paralogPairs['synthetic_fdr01'] = None;
    for index, row in paralogPairs.iterrows():
        gene1 = row[0]; gene2 = row[1];
        if gene1 in geneEss_dict:
            if geneEss_dict[gene1] == 0:
                paralogPairs.loc[index, 'synthetic'] = 0;
            elif geneEss_dict[gene1] == 1 and gene1 in dependData.columns and gene2 in genelossData.columns:
                dat1 = pd.merge(dependData[gene1],genelossData[gene2],left_index=True,right_index=True, how='inner')
                dat2 = pd.merge(dat1, cellAnnot,left_index=True,right_index=True, how='inner')
                dat2.columns = ['gene1','gene2','celltype']
                if len(dat2.gene2.unique()) > 1:
                    dat2['gene1'] = dat2['gene1'].astype(np.float64);dat2['gene2'] = dat2['gene2'].astype(np.float64);
                    ols = smf.ols('gene1 ~ C(gene2) + C(celltype)', data=dat2).fit()
                    paralogPairs.loc[index, 'coef'] = ols.params['C(gene2)[T.1.0]']
                    paralogPairs.loc[index, 'pval'] = ols.pvalues['C(gene2)[T.1.0]']
        if index % 1000 == 0:
            print("processing: ", index)

    # FDR adjust for synthetics paralog pairs
    paralogFdrAdj = paralogPairs[(paralogPairs['coef'].notna()) & (paralogPairs['pval'].notna())]
    paralogFdrAdj['FDR'] = fdrcorrection(paralogFdrAdj.pval)[1]
    paralogFdrAdj['synthetic_p005'] = np.where((paralogFdrAdj.pval < 0.05) & (paralogFdrAdj.coef < 0), 1, 0)
    paralogFdrAdj['synthetic_fdr01'] = np.where((paralogFdrAdj.FDR < 0.1) & (paralogFdrAdj.coef < 0), 1, 0)
    paralogOther = paralogPairs[~((paralogPairs['coef'].notna()) & (paralogPairs['pval'].notna()))]
    paralogPairs2 = pd.concat([paralogOther[['gene1','gene2','synthetic_p005','synthetic_fdr01']],
                               paralogFdrAdj[['gene1','gene2','synthetic_p005','synthetic_fdr01']]]) # get synthetic information ready for finalize
    
    para1df = pd.merge(para1, paralogPairs2,  how='left', left_on=['gene1','gene2'], right_on = ['gene1','gene2'])
    para2df = pd.merge(para2, paralogPairs2,  how='left', left_on=['gene1','gene2'], right_on = ['gene1','gene2'])
    paraDF = pd.merge(para1df, para2df,  how='left', left_on=['gene1','gene2'], right_on = ['gene2','gene1'])
    paraDF['synthetic_fdr_final'] = paraDF[['synthetic_fdr01_x','synthetic_fdr01_y']].max(axis=1)
    paraDF['synthetic_pval_final'] = paraDF[['synthetic_p005_x','synthetic_p005_y']].max(axis=1)
    # sort data based origanl paralog gene order
    paraDF.set_index(paraDF['gene1_x'].astype(str) + '_' + paraDF['gene2_x'].astype(str), inplace=True)
    syn_fdr_dict = paraDF['synthetic_fdr_final'].to_dict()
    syn_pval_dict = paraDF['synthetic_pval_final'].to_dict()
    input_orignial = (para1['gene1'].astype(str) + '_' + para1['gene2'].astype(str)).to_list()
    # output synthetic based both pvalue and fdr
    syn_fdr_list = [syn_fdr_dict[x] for x in input_orignial]
    syn_pval_list = [syn_pval_dict[x] for x in input_orignial]
    
    return [syn_pval_list, syn_fdr_list]


# CRISPRGeneEffect parsing
CRISPRGeneDependency = pd.read_csv(CRISPRGeneEffect,header=0,index_col=0)
CRISPRGeneDependency = CRISPRGeneDependency[CRISPRGeneDependency.columns.intersection(DepMap_gene_IDmap['ID'])].dropna(axis=1, how='any') # drop column with NaN, esp for ScreenGeneEff data
CRISPRGeneDependency = CRISPRGeneDependency.rename(columns=DepMap_gene_IDmap.set_index('ID')['Geneid']) # replace ID to ensemblid
CRISPRGeneDependency = CRISPRGeneDependency[CRISPRGeneDependency.columns.intersection(list(pd.unique(paralogData.values.ravel('K'))))] # get all paralog genes only
# only 27 cantype match with TCGA code

# paralog_synthetic(paralogPair, dependData, genelossData, geneEss_dict, cellAnnot):
Result_Synthetic = paralogData.copy();
celltypeDat = cellAnnot[['OncotreeLineage']]; celltypeDat.columns = ['celltype']
outAll = paralog_synthetic(paralogData, CRISPRGeneDependency,lossExpCnvSnv2,
                           dict(zip(geneEssential['Geneid'],geneEssential['cell_all_essential'])),
                           celltypeDat)
Result_Synthetic['cell_all_synthetic_pval'] = outAll[0]
Result_Synthetic['cell_all_synthetic_fdr'] = outAll[1]

celltypeDat = cell_reAnnot[['TCGAcode']]; celltypeDat.columns = ['celltype']
outTcga = paralog_synthetic(paralogData, CRISPRGeneDependency,lossExpCnvSnv2,
                           dict(zip(geneEssential['Geneid'],geneEssential['cell_tcga_essential'])),
                           celltypeDat)
Result_Synthetic['cell_tcga_synthetic_pval'] = outTcga[0]
Result_Synthetic['cell_tcga_synthetic_fdr'] = outTcga[1]

# Go to find essentials for per Cancer Type
for cancerType in ['BLCA','ESCA','GBM','HNSC','RCC','LIHC','NSCLC','SKCM','STAD','BRCA','COADREAD']:
    print("processing gene pair synthetic for", cancerType)
    celltypeDat0 = cell_reAnnot[cell_reAnnot['TCGAcode'] == cancerType]
    celltypeDat = celltypeDat0[['TCGAcode']]; celltypeDat.columns = ['celltype']
    out = paralog_synthetic(paralogData, CRISPRGeneDependency,lossExpCnvSnv2,
                               dict(zip(geneEssential['Geneid'],geneEssential['cell_' + cancerType + '_essential'])),
                               celltypeDat)
    Result_Synthetic['cell_' + cancerType + '_synthetic_pval'] = out[0]
    Result_Synthetic['cell_' + cancerType + '_synthetic_fdr'] = out[1]


Result_Synthetic.to_csv(str(Path(depmap_datadir)) + '/DepMap_paralogGenes_synthetic_data.csv', index=False)
#

# Func-3: paralog genes CCLE coexpression analysis
OmicsExpr = str(Path(depmap_datadir)) + "/" + version.replace(" ","-") + "_OmicsExpressionProteinCodingGenesTPMLogp1.csv"

# collect gene loss from gene expression data
exprDat = pd.read_csv(OmicsExpr,header=0,index_col=0)
exprDat = exprDat[exprDat.columns.intersection(DepMap_gene_IDmap['ID'])]
exprDat = exprDat.rename(columns=DepMap_gene_IDmap.set_index('ID')['Geneid']) # replace ID to ensemblid

def paralog_coexpression(paralogPair, geneExpData, cellAnnot):
    paralogPair.columns = ['gene1','gene2']
    rho_list = []; max_mean_vals =[]; min_mean_vals =[]; paralog_mean_vals = [];
    for index, row in paralogPair.iterrows():
        rho = None; max_mean = None; min_mean = None; paralog_mean =None;
        gene1 = row[0]; gene2 = row[1];
        if gene1 in geneExpData.columns and gene2 in geneExpData.columns:
            dat1 = geneExpData[[gene1,gene2]]
            dat2 = pd.merge(dat1, cellAnnot,left_index=True,right_index=True, how='inner')
            dat2.columns = ['gene1','gene2','celltype']
            # spearman-test
            (rho, pval) = spearmanr(dat2['gene1'], dat2['gene2'])
            
            max_mean = max(dat2['gene1'].replace('NaN', np.nan).mean(skipna=True),
                           dat2['gene2'].replace('NaN', np.nan).mean(skipna=True));
            min_mean = min(dat2['gene1'].replace('NaN', np.nan).mean(skipna=True),
                           dat2['gene2'].replace('NaN', np.nan).mean(skipna=True));
            paralog_mean = np.mean([min_mean,max_mean]);
        
        rho_list.append(rho);
        max_mean_vals.append(max_mean);
        min_mean_vals.append(min_mean);
        paralog_mean_vals.append(paralog_mean);

        if index % 1000 == 0:
            print("processing: ", index)
    return rho_list, max_mean_vals, min_mean_vals, paralog_mean_vals;

# paralog_synthetic(paralogPair, dependData, genelossData, geneEss_dict, cellAnnot):
Result_Coexpression = paralogData.copy();
celltypeDat = cellAnnot[['OncotreeLineage']]; celltypeDat.columns = ['celltype']
Result_Coexpression['cell_all_spearman.rho'],Result_Coexpression['cell_all_max_mean'],Result_Coexpression['cell_all_min_mean'],Result_Coexpression['cell_all_paralog_mean'] = paralog_coexpression(paralogData, exprDat, celltypeDat)

celltypeDat = cell_reAnnot[['TCGAcode']]; celltypeDat.columns = ['celltype']
Result_Coexpression['cell_tcga_spearman.rho'],Result_Coexpression['cell_tcga_max_mean'],Result_Coexpression['cell_tcga_min_mean'],Result_Coexpression['cell_tcga_paralog_mean']  = paralog_coexpression(paralogData, exprDat, celltypeDat)

for cancerType in ['BLCA','ESCA','GBM','HNSC','RCC','LIHC','NSCLC','SKCM','STAD', 'BRCA','COADREAD']:
    print("processing gene pair synthetic for", cancerType)
    celltypeDat0 = cell_reAnnot[cell_reAnnot['TCGAcode'] == cancerType]
    celltypeDat = celltypeDat0[['TCGAcode']]; celltypeDat.columns = ['celltype']
    Result_Coexpression['cell_' + cancerType + '_spearman.rho'],Result_Coexpression['cell_' + cancerType + '_max_mean'],Result_Coexpression['cell_' + cancerType + '_min_mean'],Result_Coexpression['cell_' + cancerType + '_paralog_mean'] = paralog_coexpression(paralogData, exprDat, celltypeDat)

Result_Coexpression.to_csv(str(Path(depmap_datadir)) + '/Paralog_DepMap_co-expression_data.csv', index=False)
#
