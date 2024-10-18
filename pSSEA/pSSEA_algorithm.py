#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 2023

@author: cpdong
"""
import argparse;
import numpy as np
import pandas as pd
import time;

parser = argparse.ArgumentParser()
parser.add_argument('-input', metavar = 'input data', dest='input', help='input file');
parser.add_argument('-output', metavar = 'output data', dest='output', help='output file');
args = parser.parse_args();


"""
sg_pergene: sgrna per gene, typically 4 -5, if variable, pick the larget one
nperm: number of permulation for calculating p value
"""

def cdf(sample, x, sort = False):
    # Sorts the sample, if unsorted
    if sort:
        sample.sort()
    # Counts how many observations are below x
    cdf = sum(sample <= x)
    # Divides by the total number of observations
    cdf = cdf / len(sample)
    return cdf

def ks_test(gset, geneRanking):
    """
    gset: genes belongs to a gene sets
    geneRanking: ranked all sgrna_ids based lfc or p-value
    """
    tag_indicator = np.in1d(geneRanking, gset, assume_unique=True).astype(int)
    gset_idxs = np.flatnonzero(tag_indicator)
    offset_idxs = np.where(tag_indicator == 0)[0]
    
    # Evaluates the KS statistic
    statistics = [] # KS Statistic list
    for x in gset_idxs:
        cdf_gset = cdf(sample = gset_idxs, x  = x)
        cdf_offset = cdf(sample = offset_idxs, x  = x)
        statistics.append(cdf_gset - cdf_offset)
        
    # https://towardsdatascience.com/decoding-gene-set-variation-analysis-8193a0cfda3
    ks_stat = max(0, max(statistics)) - abs(min(0, min(statistics))) 
    return ks_stat


def ks_test_null(NTC_sets, gsetSize, geneRanking, nperm = 10000):
    """
    pre-compute all null rndwalk stat
    NTC_sets:  all NTC sgrnas
    gsetSize: genes belongs to a gene sets
    geneRanking: ranked all sgrna_ids based lfc or p-value or lfc
    """
    
    tag_indicator = np.concatenate([np.full(gsetSize, 1), np.full(len(NTC_sets) - gsetSize, 0)]); # create initial array
    tag_indicators = np.tile(tag_indicator, (nperm, 1))
    
    rs = np.random.RandomState(123)
    for i in range(nperm):
        rs.shuffle(tag_indicators[i])
    
    ks_stat_null = []
    for i in range(len(tag_indicators)):
        ntc_set = list(NTC_sets[np.flatnonzero(tag_indicators[i])]); # get random NTC sgrnas
        tag_indicator = np.in1d(geneRanking, ntc_set, assume_unique=True).astype(int)
        gset_idxs = np.flatnonzero(tag_indicator)
        offset_idxs = np.where(tag_indicator == 0)[0]
        
        statistics = [] # KS statistic list
        for x in gset_idxs:
            cdf_gset = cdf(sample = gset_idxs, x  = x)
            cdf_offset = cdf(sample = offset_idxs, x  = x)
            statistics.append(cdf_gset - cdf_offset)
            
        ks_stat_null.append(max(0, max(statistics)) - abs(min(0, min(statistics))))
    return np.array(ks_stat_null)

def ks_test_null_arrays(gsetSize_min, gsetSize_max, NTC_sets, geneRanking, nperm = 10000): # precompute all null between min - max gset size
    ks_null_list = [None] * (gsetSize_max + 1); # pre-build array for storaging
    for i in range(gsetSize_max + 1):
        if i >= gsetSize_min and i <= gsetSize_max:
            # ks_test_null(NTC_sets, gsetSize, geneRanking, nperm = 10000)
            ks_null_list[i] = ks_test_null(NTC_sets, i+1, geneRanking, nperm = 10000);
        elif i > gsetSize_max:
            break;
            
    return ks_null_list

def compute_pval(es, esnull):
    """
    Adopt from https://github.com/zqfang/GSEApy/tree/master/gseapy
    output pval: normimal pvalue
    """
    condlist = [es < 0, es >= 0]
    choicelist = [
        len(esnull[esnull < es]) / max(len(esnull[esnull < 0]), len(esnull[esnull < es])),
        len(esnull[esnull > es]) / max(len(esnull[esnull > 0]), len(esnull[esnull > es]))
    ]
    pvals = np.select(condlist, choicelist)

    return pvals

if __name__ == "__main__":
    
    """
    General K-S test settings
    """
    nperm = 10000
    Set_Max = 10;
    Set_Min = 4;
    
    
    #DEGfile = '/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/10_label/GSE149933_data/GSE149933_DESeq2_TKO_B16_result.tsv'
    DEGfile = args.input;
    result = pd.read_csv(DEGfile, header=0, sep='\t')
    result['symbol'] = result['sgrna'].str.split('_', expand=True).iloc[:, 1]
    
    NTC_sets = np.array(result[result['NTC']==1]['tmpid']); # get NTC sgrna list

    annot = pd.read_csv("https://raw.githubusercontent.com/cpdong/ParalogICB/main/data/MouseHuman_matched_annotation.tsv", header=0, sep="\t")
    annot = annot[['ensembl_gene_id','hgnc_symbol','mgi_symbol']]
    
    paralogData = pd.read_csv('/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/10_label/00_Ensembl102_pcgene_paralog_bioMart_query_result.tsv', header=0, sep='\t')
    paralogData['newID'] = [ '-'.join(sorted(x)) for x in paralogData[["ensembl_gene_id","hsapiens_paralog_ensembl_gene"]].values.tolist() ]
    paralogData = paralogData.drop_duplicates(subset='newID', keep="first")
    paralogData = paralogData[["ensembl_gene_id","hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')
    
    result2 = pd.merge(annot, result, left_on='mgi_symbol', right_on='symbol', how='inner')
    result2 = result2[result2['ensembl_gene_id'].isin(list(pd.unique(paralogData.values.ravel('K'))))]; # get all gene in the paralogData
    freq = result2["ensembl_gene_id"].value_counts()
    filtered = freq[freq<= 5].index.values.tolist() # only pick gene with sgRNAs <=5
    result3 = result2[result2['ensembl_gene_id'].isin(filtered)];
    result3 = result3.drop_duplicates(subset='tmpid', keep="first")
    #result3.to_csv('/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/10_label/GSE149933_data/GSE149933_DESeq2_TKO_B16_humanized_data.tsv', index=False, sep='\t');
    
    # get NTC data and rbind with select ones
    result_NTC = result[result['NTC']==1];
    result_NTC['ensembl_gene_id'] = result_NTC['hgnc_symbol'] = result_NTC['mgi_symbol'] = result_NTC['symbol'] ='NTC';

    result4 = pd.concat([result3, result_NTC])
    result5 = result4.sort_values(by=['log2FoldChange'], ascending=True) # depletion test with low to high
    
    geneRanking = result5['tmpid'].values.tolist()
    
    '''
    runing ES with KS test
    pre-compute null test
    '''
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'KS-test start!')
    ks_test_null_list = ks_test_null_arrays(gsetSize_min=4,
                                              gsetSize_max=10, 
                                              NTC_sets=NTC_sets,
                                              geneRanking=geneRanking,
                                              nperm = 10000)

    sgrnas_total = []; sgrnas_g1 = []; sgrnas_g2 = [];
    ES_ks = []; ES_ks_pval = [];
    
    count_per_gene = result5['ensembl_gene_id'].value_counts().to_dict()
    sgrna_per_gene = result5.groupby('ensembl_gene_id')['tmpid'].apply(list).to_dict()
    for index, row in paralogData.iterrows():
        gene1 = row[0]; gene2 = row[1];
        sgrna_g1 = sgrna_g2 = sgrna_total = rndwalk_stat = rndwalk_pval = ks_stat = ks_pval = None;
        if gene1 in result5['ensembl_gene_id'].values.tolist() and gene2 in result5['ensembl_gene_id'].values.tolist():
            sgrna_g1 = count_per_gene[gene1]
            sgrna_g2 = count_per_gene[gene2]
            sgrna_set = sgrna_per_gene[gene1] + sgrna_per_gene[gene2]
            sgrna_total = sgrna_g1 + sgrna_g2;
            if sgrna_total >= Set_Min and sgrna_total <= Set_Max:
                ks_stat = ks_test(gset=sgrna_set, geneRanking=geneRanking)
                ks_pval = compute_pval(ks_stat, ks_test_null_list[sgrna_total])
                
                
        sgrnas_g1.append(sgrna_g1)
        sgrnas_g2.append(sgrna_g2)
        sgrnas_total.append(sgrna_total)
        ES_ks.append(ks_stat)
        ES_ks_pval.append(ks_pval)
        
    paralogData['sgrnas_g1'] = sgrnas_g1
    paralogData['sgrnas_g2'] = sgrnas_g2
    paralogData['sgrnas_total'] = sgrnas_total
    paralogData['ES_ks'] = ES_ks
    paralogData['ES_ks_pval'] = ES_ks_pval
    
    paralogData.to_csv(args.output, index=False, sep='\t');
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'KS-test End!')
    #
