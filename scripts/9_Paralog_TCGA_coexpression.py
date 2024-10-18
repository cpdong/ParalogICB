import gzip, os, shutil, time;
import pandas as pd;
from pathlib import Path
from scipy.stats import spearmanr
import urllib.request;


'''
Use brew x86 to install x86_64 python 3.9  on M1 Chip MacOS, ref https://www.qiniu.com/qfans/qnso-70315418#comments
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
alias brew86="arch -x86_64 /usr/local/bin/brew" # add to zshrc file
brew86 install python@3.9
'''

# synthetic gene pairs for survival: 10.1016/j.celrep.2019.06.067
# new scores: Mapping the landscape of synthetic lethal interactions in liver cancer: Functional Similarity Analysis MF and CC


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

# synthetic gene pairs for survival: 0.1016/j.celrep.2019.06.067
# new scores: Mapping the landscape of synthetic lethal interactions in liver cancer

if __name__ == "__main__":
    
    workdir="/Users/cpdong/Downloads/tcga/"
    
    paralog_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')
    
    
    # Go through ALL,NSCLC,RCC,BLCA,ESCA,GBM,HNSC,LIHC,SKCM,STAD
    coexpr_TCGA  = paralogData.copy()
    cantypes = ['ALL', 'RCC','SKCM','BLCA','STAD','LIHC','HNSC','ESCA','GBM','NSCLC']
    for cancer in cantypes:
        cancer = 'GBM'
        if cancer == 'ALL':
            tcgaFile = str(Path(workdir)) + '/TCGA-PANCAN_tumor_tpm.tsv' # note the TCGA file path
        else:
            tcgaFile = str(Path(workdir)) + '/TCGA-' + cancer + '_tumor_tpm.tsv' # note the TCGA file path

        tcgaData = pd.read_csv(tcgaFile, header=0, index_col=0, sep="\t").T
            
        rhos = []
        for index, row in coexpr_TCGA.iterrows():
            rho = None; gene1 = row[0]; gene2 = row[1];
            if gene1 in tcgaData.columns and gene2 in tcgaData.columns:
                (rho, pval) = spearmanr(tcgaData[gene1], tcgaData[gene2]) # coexpression function
            
            if index % 1000 == 0:
                print("processing: ", cancer, index)
            rhos.append(rho)
        
        coexpr_TCGA[cancer] = rhos
    
    coexpr_TCGA.to_csv(str(Path(workdir)) + '/paralogGenes_co-expression_TCGA_data.csv', index=False)
    del tcgaData, coexpr_TCGA;
    
    
