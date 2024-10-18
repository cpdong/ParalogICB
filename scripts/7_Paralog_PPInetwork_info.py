import argparse, os, re;
import gzip,zipfile,shutil, time;
import bs4,subprocess;
import requests
import pandas as pd;
import numpy as np;
from pathlib import Path
from functools import reduce
import scipy.stats as stats
from scipy.stats import kendalltau
from scipy.stats import spearmanr
from scipy.stats import zscore
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.formula.api as smf
import urllib.request;


def BioGRID_download(savedir='./'):
    if savedir=='./':
        savedir = os.getcwd()
    
    #savedir = '/Users/cpdong/Downloads/'
    if not os.path.exists(str(Path(savedir)) + '/BioGRID'):
        os.makedirs(str(Path(savedir)) + '/BioGRID')
    # database
    release_url = 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.225/BIOGRID-ORGANISM-4.4.225.tab.zip'
    urllib.request.urlretrieve(release_url, str(Path(savedir)) + '/BioGRID/BIOGRID-ORGANISM-4.4.225.tab.zip');
    with zipfile.ZipFile(str(Path(savedir)) + '/BioGRID/BIOGRID-ORGANISM-4.4.225.tab.zip', 'r') as zipf:
        zipf.extractall(str(Path(savedir)) + '/BioGRID')
    shutil.copy2(str(Path(savedir)) + '/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-4.4.225.tab.txt', savedir)
    
    # database identifiers and processing
    release_id_url =  'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.225/BIOGRID-IDENTIFIERS-4.4.225.tab.zip'
    urllib.request.urlretrieve(release_id_url, str(Path(savedir)) + '/BioGRID/BIOGRID-IDENTIFIERS-4.4.225.tab.zip');
    with zipfile.ZipFile(str(Path(savedir)) + '/BioGRID/BIOGRID-IDENTIFIERS-4.4.225.tab.zip', 'r') as zipf:
        zipf.extractall(str(Path(savedir)) + '/BioGRID')
    
    with open(str(Path(savedir)) + '/BioGRID/BIOGRID-IDENTIFIERS-4.4.225.tab.txt') as f1:
        with open(str(Path(savedir)) + '/BioGRID/BIOGRID-IDENTIFIERS.txt', 'w') as f2:
            lines = f1.readlines()
            for line in lines:
                if 'Homo sapiens' in line:
                    f2.write(line)
                    
    d = pd.read_csv(str(Path(savedir)) + '/BioGRID/BIOGRID-IDENTIFIERS.txt', header=None, sep='\t')
    d1 = d.loc[d.iloc[:,2] =="OFFICIAL SYMBOL", 0:1];d1.columns = ['entrez_gene_id','hgnc_symbol']
    d2 = d.loc[d.iloc[:,2] =="ENSEMBL", 0:1]; d2.columns = ['entrez_gene_id','ensemble_gene_id']
    d3 = pd.merge(d1, d2, left_on='entrez_gene_id', right_on='entrez_gene_id', how='inner')
    d3.to_csv(str(Path(savedir)) + '/BIOGRID-IDs.txt', index=False, sep='\t')
    shutil.rmtree(str(Path(savedir)) + '/BioGRID') # remnove the temp folder


def compute_ppi_summary_for_pairs(paralogPairs, ppi_data, ppi_idmap, geneEss):
    # Make ppi df symmetrical for merging with paralog pairs (which are the unique list)
    # ppi_data = ppi_data[ppi_data['Experimental System Type']=='physical']; # if you want pick physical link
    ppi_data = BioGRID_db.copy(); ppi_idmap = BioGRID_idmap.copy(); paralogPairs = paralogData.copy();
    paralogPairs.columns = ['gene1','gene2']
    p1 = pd.merge(ppi_data, ppi_idmap, left_on='OFFICIAL_SYMBOL_A', right_on='hgnc_symbol', how ='inner')
    p2 = pd.merge(p1, ppi_idmap, left_on='OFFICIAL_SYMBOL_B', right_on='hgnc_symbol', how ='inner')
    
    def infer_gene_interact(paralogPairs, ppi_data): # gene - gene interaction inferrings
        # process BioGRID data
        # ppi_data = ppi_data[ppi_data['Experimental System Type']=='physical']; # only pick physical link
        v1 = ppi_data[['ensemble_gene_id_x', 'ensemble_gene_id_y']].values.tolist();
        ppi_data['psuedoid'] = ['-'.join(sorted(x)) for x in v1]; #sort element
        # process paralogs data
        paralogPairs['psuedoid'] =  [ '-'.join(sorted(x)) for x in paralogPairs.values.tolist() ]
        
        return list(np.where(paralogPairs['psuedoid'].isin(p2['psuedoid']), 1, 0))

    paralogPairs['interact'] = infer_gene_interact(paralogPairs= paralogPairs, ppi_data= p2); # get gene interact result
    
    # Following function adapt from: https://github.com/cancergenetics/paralog_SL_prediction
    ppi_symmetric = pd.concat([p2[['ensemble_gene_id_x', 'ensemble_gene_id_y']],
                               p2[['ensemble_gene_id_y', 'ensemble_gene_id_x']]]).drop_duplicates(keep='first')
    ppi_symmetric = ppi_symmetric.reset_index(drop=True)
    
    # Use symmetric version of ppi table to get set of all interactors for each gene
    ppi_per_gene = ppi_symmetric.groupby('ensemble_gene_id_x').agg({'ensemble_gene_id_y':set}).reset_index().rename(columns={'ensemble_gene_id_x':'gene', 'ensemble_gene_id_y':'ppi'})

    # Merge ppi_per_gene with each A1 and A2 in all pairs
    # Note: pairs can have shared interactors even if there is no evidence they interact themselves
    df = pd.merge(paralogPairs, ppi_per_gene.rename(columns={'gene':'gene1','ppi':'A1_ppi'}), on='gene1', how='left')
    df = pd.merge(df, ppi_per_gene.rename(columns={'gene':'gene2','ppi':'A2_ppi'}), on='gene2', how='left')
    
    # Fill NaNs with empty sets
    df['A1_ppi'] = df['A1_ppi'].apply(lambda d: d if not pd.isnull(d) else set())
    df['A2_ppi'] = df['A2_ppi'].apply(lambda d: d if not pd.isnull(d) else set())
    
    # Remove A2 gene in the set of interactors for A1 gene (and vice versa)
    # Don't want to include these in union for other calculations
    df.A1_ppi = df.apply(lambda x: x.A1_ppi - {x.gene2}, axis=1)
    df.A2_ppi = df.apply(lambda x: x.A2_ppi - {x.gene1}, axis=1)

    # Calculate total num interactors + shared interactors
    df['n_A1_ppi'] = df.apply(lambda x: len(x.A1_ppi), axis=1)
    df['n_A2_ppi'] = df.apply(lambda x: len(x.A2_ppi), axis=1)
    df['shared_ppi'] = df.apply(lambda x: x.A1_ppi.intersection(x.A2_ppi), axis=1)
    df['n_total_ppi'] = df.apply(lambda x: len(x.A1_ppi.union(x.A2_ppi)), axis=1)
    df['n_shared_ppi'] = df.apply(lambda x: len(x.shared_ppi), axis=1)                

    # Calculate jaccard index for shared interactors
    def calc_jaccard_index(x):
        if x.n_shared_ppi == 0: return 0
        return x.n_shared_ppi / ((x.n_A1_ppi + x.n_A2_ppi) - x.n_shared_ppi)

    df['shared_ppi_jaccard_idx'] = df.apply(calc_jaccard_index, axis=1)

    # Calculate FET for overlap of interactors, N = all genes involved in interactions
    N = len(pd.concat([paralogPairs.gene1, paralogPairs.gene2]).unique())
    print('N genes involded in paralog pairs data:', N)
    # ctab:   | A2      | Not A2
    #      A1 | shared  | A1 only
    #  Not A1 | A2 only | N - union(A1, A2)
    def calc_fet_shared_ppi(x):
        ctab = pd.DataFrame({'A2': [x.n_shared_ppi, x.n_A2_ppi - x.n_shared_ppi],
                            'NA2': [x.n_A1_ppi - x.n_shared_ppi, N - x.n_total_ppi]}, index=['A1', 'NA1'])
        (OR, pval) = stats.fisher_exact(ctab)
        if pval==0: # Use smallest float64 number, to apply log10
            pval = np.nextafter(0, 1)
        log_pval = (-np.log10(pval)) if pval != 1 else 0
        log_pval = -log_pval if OR < 1 else log_pval # return negative pval if signif depletion in overlap
        return log_pval

    df['fet_ppi_overlap'] = df.apply(calc_fet_shared_ppi, axis=1)
    
    # Calculate essential of shared interactors
    Essential_perc = []; Essential_Chronos = [];
    for index, row in df.iterrows():
        Ess_perc = None; Ess_chronos_avg = None;
        shared_ppi = row['shared_ppi'];
        if len(shared_ppi) > 0:
            geneEss_tmp = geneEss[geneEss.index.isin(shared_ppi)]; 
            
            Ess_perc = geneEss_tmp['cell_all_essential'].replace('NaN', np.nan).mean(skipna=True);
            Ess_chronos_avg = geneEss_tmp['cell_all_avgChronos'].replace('NaN', np.nan).mean(skipna=True);
    
        Essential_perc.append(Ess_perc);
        Essential_Chronos.append(Ess_chronos_avg);
    
    df['shared_ppi_Ess_perc'] = Essential_perc;
    df['shared_ppi_Ess_chronos_avg'] = Essential_Chronos;
    df = df.drop(columns=['psuedoid','A1_ppi','A2_ppi','n_A1_ppi','n_A2_ppi','shared_ppi'])    
    return df


def compute_colocalize(paralogPairs, tmpdir='./'):
    if tmpdir=='./':
        tmpdir = os.getcwd()
    
    # tmpdir='/Users/cpdong/Downloads/'
    colocal_url = 'https://www.proteinatlas.org/download/subcellular_location.tsv.zip'
    urllib.request.urlretrieve(colocal_url, str(Path(tmpdir)) + '/colocalized.zip');
    with zipfile.ZipFile(str(Path(tmpdir)) + '/colocalized.zip', 'r') as zipf:
        zipf.extractall(Path(tmpdir))
    df = pd.read_csv(str(Path(tmpdir)) + '/subcellular_location.tsv',header=0,sep='\t')
    df['locations'] = df['Main location'] + ';' + df['Additional location']
    df['locations'] = df['locations'].apply(lambda d: set(d.split(';')) if not pd.isnull(d) else set())
    df = df[['Gene','locations']]
    
    # paralogPairs = paralogData.copy();
    paralogPairs.columns = ['gene1','gene2']
    pdat = pd.merge(paralogPairs, df.rename(columns={'Gene':'gene1','locations':'A1_location'}), on='gene1', how='left')
    pdat = pd.merge(pdat, df.rename(columns={'Gene':'gene2','locations':'A2_location'}), on='gene2', how='left')
    
    # Calculate total num interactors + shared interactors
    pdat['A1_location'] = pdat['A1_location'].apply(lambda d: d if not pd.isnull(d) else set())
    pdat['A2_location'] = pdat['A2_location'].apply(lambda d: d if not pd.isnull(d) else set())
    pdat['n_shared_locals'] = pdat.apply(lambda x: len(x.A1_location.intersection(x.A2_location)), axis=1)
    pdat['n_total_locals'] = pdat.apply(lambda x: len(set(x.A1_location.union(x.A2_location))), axis=1)
    colocalize_jaccard_idx = pdat.apply(lambda x: x.n_shared_locals / x.n_total_locals if x.n_total_locals !=0 else None, axis=1)
    os.remove(str(Path(tmpdir)) + '/colocalized.zip') # remnove the temp file
    os.remove(str(Path(tmpdir)) + '/subcellular_location.tsv') # remnove the temp file
    
    return colocalize_jaccard_idx

def complex_CORUM(paralogPairs, geneEss, tmpdir='./'):
    if tmpdir=='./':
        tmpdir = os.getcwd()
    
    tmpdir='/Users/cpdong/Downloads/'
    CORUM_url = 'https://maayanlab.cloud/static/hdfs/harmonizome/data/corum/gene_attribute_matrix.txt.gz'
    urllib.request.urlretrieve(CORUM_url, str(Path(tmpdir)) + '/CORUM_db.gz');
    with gzip.open(str(Path(tmpdir)) + '/CORUM_db.gz', 'r') as fin, open(str(Path(tmpdir)) + '/CORUM_db', 'wb') as fout:
        shutil.copyfileobj(fin, fout)
    df = pd.read_csv(str(Path(tmpdir)) + '/CORUM_db',header=0,skiprows=[1,2], sep='\t')
    df.drop(df.columns[[0,1]], axis=1, inplace=True)
    
    annot = pd.read_csv('https://raw.githubusercontent.com/cpdong/ParalogICB/main/data/gencode.v36_annotation.tsv',header=0,sep='\t')
    annot = annot.loc[annot['entrezgene_id'].notna(), ['Geneid', 'entrezgene_id']]
    annot['entrezgene_id'] = annot['entrezgene_id'].astype(int)
    
    df = pd.merge(annot, df, left_on='entrezgene_id', right_on='Complex', how='inner')
    df = df.drop(columns=['entrezgene_id', 'Complex'])
    df.set_index('Geneid', inplace = True)
    df2 = df.T
    
    both_in_same_complex = []; both_in_complex = []; Essential_perc = []; Essential_Chronos = [];# but no need to be in the same one complex
    for index, row in paralogPairs.iterrows():
        in_complex = 0; in_same_complex = 0; Ess_perc = None; Ess_chronos_avg = None;
        gene1 = row[0]; gene2 = row[1];
        if any(x in df2.columns for x in [gene1, gene2]):
            tmpdf =  df2[df2.columns.intersection([gene1,gene2])]
            tmpdf['sums'] = tmpdf.sum(axis=1, numeric_only=True)
            # tmpdf['sums'] = df[gene1] + df[gene2];
            if tmpdf['sums'].max() >= 1:
                in_complex = 1;
                if tmpdf['sums'].max() == 2:
                    in_same_complex = 1;
                
                complexNames=tmpdf.index[tmpdf['sums'] >= 1].tolist();
                complDat = df[df.columns.intersection(complexNames)];
                complDat['sums'] = complDat.sum(axis=1, numeric_only=True);
                complex_all_genes = complDat.index[complDat['sums'] >= 1].tolist();
                geneEss_tmp = geneEss[geneEss.index.isin(complex_all_genes)];
                
                Ess_perc = geneEss_tmp['cell_all_essential'].replace('NaN', np.nan).mean(skipna=True);
                Ess_chronos_avg = geneEss_tmp['cell_all_avgChronos'].replace('NaN', np.nan).mean(skipna=True);
                
        both_in_same_complex.append(in_same_complex);
        both_in_complex.append(in_complex);
        Essential_perc.append(Ess_perc);
        Essential_Chronos.append(Ess_chronos_avg);
        
    os.remove(str(Path(tmpdir)) + '/CORUM_db.gz') # remnove the temp file
    os.remove(str(Path(tmpdir)) + '/CORUM_db') # remnove the temp file 
    return both_in_same_complex, both_in_complex, Essential_perc, Essential_Chronos;


if __name__ == "__main__":
    
    workdir="/path/to/workdir"
    paralog_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    paralogData = pd.read_csv(paralog_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id", "hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')

    tmpdir = './tmp/'
    BioGRID_download(savedir = tmpdir); # download biogrid data and ids
    BioGRID_db = pd.read_csv(str(Path(tmpdir)) + '/BIOGRID-ORGANISM-Homo_sapiens-4.4.225.tab.txt', header=0, sep='\t', skiprows=35)
    BioGRID_idmap = pd.read_csv(str(Path(tmpdir)) + '/BIOGRID-IDs.txt', header=0, sep='\t')[['hgnc_symbol','ensemble_gene_id']]

    geneEss_data = pd.read_csv('./DepMap_data/DepMap_singleGene_essential_avgChronos_data.csv', header=0, index_col=0)[['cell_all_essential','cell_all_avgChronos']]
    df = compute_ppi_summary_for_pairs(paralogPairs=paralogData, ppi_data=BioGRID_db, ppi_idmap=BioGRID_idmap, geneEss=geneEss_data)
    
    df['colocalization'] = compute_colocalize(paralogPairs=paralogData.copy(), tmpdir=tmpdir)
    df['both_in_same_complex'], df['in_CORUM_complex'], df['complex_Ess_perc'], df['complex_Ess_Chronos'] = complex_CORUM(paralogPairs=paralogData.copy(), geneEss=geneEss_data, tmpdir='/Users/cpdong/Downloads/')
    
    df.to_csv(str(Path(workdir)) + '/Paralog_ppi_colocalize_complex_information.tsv', index=False, sep='\t')
    
    os.remove(str(Path(tmpdir)) + '/BIOGRID-IDs.txt')
    os.remove(str(Path(tmpdir)) + '/BIOGRID-ORGANISM-Homo_sapiens-4.4.225.tab.txt')
    
