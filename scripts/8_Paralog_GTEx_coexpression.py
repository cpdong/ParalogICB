import gzip, os, shutil, time;
import pandas as pd;
import numpy as np;
from pathlib import Path
from scipy.stats import spearmanr
import urllib.request;

'''
Use brew x86 to install x86_64 python 3.9  on M1 Chip MacOS, ref https://www.qiniu.com/qfans/qnso-70315418#comments
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
alias brew86="arch -x86_64 /usr/local/bin/brew" # add to zshrc file
brew86 install python@3.9
'''

# synthetic gene pairs for survival: 0.1016/j.celrep.2019.06.067
# new scores: Mapping the landscape of synthetic lethal interactions in liver cancer


def GTEx_download(savedir='./'):
    '''
    savedir: a directory for storage temp/final dowmloading data
    '''
    if savedir=='./':
        savedir = os.getcwd()
    
    #savedir = '/Users/cpdong/Downloads/tcga/tttt/'
    dl_url = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
    urllib.request.urlretrieve(dl_url, str(Path(savedir)) + '/GTEx_v8_gene_tpm.gct.gz'); # download GTEx v8 TPM data
    with gzip.open(str(Path(savedir)) + '/GTEx_v8_gene_tpm.gct.gz', 'r') as fin, open(str(Path(savedir)) + '/GTEx_v8_gene_tpm.gct', 'wb') as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(str(Path(savedir)) + '/GTEx_v8_gene_tpm.gct.gz');
    anno_url = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
    urllib.request.urlretrieve(anno_url, str(Path(savedir)) + '/GTEx_v8_Annotation.txt'); # download GTEx v8 TPM data

if __name__ == "__main__":
    
    workdir="/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/08_features/"
    os.chdir(workdir)
    
    paralog_file="/Users/cpdong/Library/CloudStorage/Dropbox/project/paralogs/data/07_coevolute//00_Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

    # download GTEx data
    if not os.path.exists('gtex'):
        os.makedirs('gtex')
    GTEx_download(savedir=str(Path(workdir)) + '/gtex')
    
    gtex_file = str(Path(workdir)) + '/gtex/GTEx_v8_gene_tpm.gct'
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ": start read!")
    gtex = pd.read_csv(gtex_file, header = 0, sep = "\t", skiprows=2)
    del gtex['Description'];
    gtex = gtex[~gtex['Name'].str.contains("_PAR")] # remove abnormal genes
    gtex.set_index('Name', inplace=True)
    gtex.index = [ x[:15] for x in gtex.index ]
    gtexData = gtex[gtex.index.isin(pd.unique(paralogData[['ensembl_gene_id','hsapiens_paralog_ensembl_gene']].values.ravel('K')))]
    del gtex
    
    GTEx_tissue_dict ={'GBM':'Brain', 'SKCM':'Skin', 'HNSC':'Salivary Gland','NSCLC':'Lung','ESCA':'Esophagus',
                       'STAD':'Stomach', 'LIHC':'Liver', 'RCC':'Kidney', 'BLCA':'Bladder', 'BRCA':'Breast', 'COADREAD':'Colon'}
    annot = pd.read_csv(str(Path(workdir)) + '/gtex/GTEx_v8_Annotation.txt', header=0, sep="\t")
    
    
    # Go through ALL,NSCLC,RCC,BLCA,ESCA,GBM,HNSC,LIHC,SKCM,STAD
    coexpr_GTEx  = paralogData.copy()
    cantypes = ['ALL', 'RCC','SKCM','BLCA','STAD','LIHC','HNSC','ESCA','GBM','NSCLC', 'BRCA', 'COADREAD']
    for tissue in cantypes:
        if tissue == 'ALL':
            gdata = gtexData.copy().T;
        else:
            tissue_samples = annot[annot['SMTS'] == GTEx_tissue_dict[tissue]]['SAMPID'].tolist()
            gdata = gtexData[gtexData.columns.intersection(tissue_samples)].T
        rho_list = []; max_mean_vals =[]; min_mean_vals =[]; paralog_mean_vals = [];
        for index, row in coexpr_GTEx.iterrows():
            rhos = None; max_mean = None; min_mean = None; paralog_mean =None;
            gene1 = row[0]; gene2 = row[1];
            if gene1 in gdata.columns and gene2 in gdata.columns:
                (rho, pval) = spearmanr(gdata[gene1], gdata[gene2]) # coexpression function
                
                max_mean = max(gdata[gene1].replace('NaN', np.nan).mean(skipna=True),
                               gdata[gene2].replace('NaN', np.nan).mean(skipna=True));
                min_mean = min(gdata[gene1].replace('NaN', np.nan).mean(skipna=True),
                               gdata[gene2].replace('NaN', np.nan).mean(skipna=True));
                paralog_mean = np.mean([min_mean,max_mean]);
                
            if index % 1000 == 0:
                print("processing: ", tissue, index)
            rho_list.append(rho)
            max_mean_vals.append(max_mean);
            min_mean_vals.append(min_mean);
            paralog_mean_vals.append(paralog_mean);
        
        coexpr_GTEx[tissue] = rho_list;
        coexpr_GTEx[tissue + '_max_mean'] = max_mean_vals;
        coexpr_GTEx[tissue + '_min_mean'] = min_mean_vals;
        coexpr_GTEx[tissue + '_paralog_mean'] = paralog_mean_vals;
    
    coexpr_GTEx.to_csv(str(Path(workdir)) + '/Paralog_GTEx_co-expression_data.csv', index=False)
    del gdata, gtexData, coexpr_GTEx;
    shutil.rmtree(str(Path(workdir)) + '/gtex') # remnove the temp folder
    
