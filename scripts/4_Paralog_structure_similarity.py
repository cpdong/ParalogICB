#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:21:25 2023

@author: cpdong
"""

import os
import re
import gzip
import time
import subprocess
import shutil
import pandas as pd
import urllib.request
import multiprocessing
import warnings
import tarfile
from pathlib import Path
from Bio import SeqIO, BiopythonWarning  # uniprotIDchecking

# ~ 24*6 hrs
outdir = "/path/to/workdir"
os.chdir(outdir)
thread = 48

def get_ensembl_uniprot_dict(pdb_dir):
    uniprot_idmapping_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
    # storage af_pdb and its prot seq length
    pdblist = os.listdir(pdb_dir)
    pdb_dict = {}
    for pfile in pdblist:
        if '-model' in pfile and ".pdb" in pfile:
            pid = re.search('AF-(.*)-F', pfile).group(1)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', BiopythonWarning)
                open_method = gzip.open if pfile.lower().endswith('.gz') else open
                with open_method(str(Path(pdb_dir)) + "/" + pfile, "rt") as pdb:
                    seqs = ''
                    for r in SeqIO.parse(pdb, "pdb-atom"):
                        seqs += str(r.seq)
            pdb_dict[pid] = len(seqs)

    # Build the dict of ensembl id - pdb abspath
    urllib.request.urlretrieve(uniprot_idmapping_url, str(Path(pdb_dir)) + "/uniprot_idmapping.dat.gz"); # download uniprot id mapping
    # open_method = gzip.open if uniprot_hsa_idmapping_dat.lower().endswith('.gz') else open;
    with gzip.open(str(Path(pdb_dir)) + "/uniprot_idmapping.dat.gz", "rt") as dat:
        dicts = {}
        for line in dat:
            if 'Ensembl' in line and 'ENSG' in line:
                accessid = line.strip().split('\t')[0]
                ensemblid = line.strip().split('\t')[2]
                if accessid in pdb_dict:
                    pdbname = [x for x in pdblist if accessid in x][0]
                    if ensemblid not in dicts:
                        dicts[ensemblid] = str(Path(pdb_dir)) + pdbname
                    else:
                        len0 = pdb_dict[dicts[ensemblid]]  # exist pdb length
                        len1 = pdb_dict[accessid]  # current pdb length
                        if len1 > len0:  # update with the longest prot
                            dicts[ensemblid] = str(Path(pdb_dir)) + pdbname

    os.remove(str(Path(pdb_dir)) + "/uniprot_idmapping.dat.gz")
    return dicts


def tmscore(gene1, gene2, prot_a, prot_b):
    TMscore = "/Users/cpdong/Documents/Biosoft/TMscore/TMscore"
    rmsd = tm_score = maxsub = 'NA'
    if os.path.isfile(prot_a) and os.path.isfile(prot_b):
        out = subprocess.Popen([TMscore, prot_a, prot_b], shell=False,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
        data = str(out).split("\\n")
        for d in data:
            x = re.sub(r"\s\s+", " ", d).split(' ')
            if x[0] == "TM-score" and x[1] == "=":
                tm_score = float(x[2])
            elif x[0] == "RMSD":
                rmsd = float(x[5])
            elif x[0] == "MaxSub-score=":
                maxsub = float(x[1])
    else:
        raise Exception("Check that %s and %s exists" % (prot_a, prot_b))
    return [gene1, gene2, rmsd, tm_score, maxsub]


if __name__ == '__main__':

    # Download alphafold predicted structure pdbs
    if not os.path.exists(str(Path(outdir)) + "/AF2pdb"):
        os.makedirs(outdir + "/AF2pdb")
    urllib.request.urlretrieve("https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar",
                               str(Path(outdir)) + "/AF2pdb.tar")
    tar = tarfile.open(str(Path(outdir)) + "/AF2pdb.tar")
    for member in tar.getmembers():
        if '.pdb' in member.name:
            tar.extract(member, path=outdir + "/AF2pdb")
            if (member.name).endswith('.gz'):  # check for ".gz" extension
                with gzip.open(outdir + "/AF2pdb/" + member.name, 'rb') as f_in, open(outdir + "/AF2pdb/" + member.name.replace('.gz', ''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                # delete gziped file
                os.remove(outdir + "/AF2pdb/" + member.name)
    tar.close()
    os.remove(str(Path(outdir)) + "/AF2pdb.tar")  # delete downloaded DB file

    ensembl_uniprot_dict = get_ensembl_uniprot_dict(outdir + "/AF2pdb/")

    ensembl102_paralog_seq_file = "Ensembl102_paralog_pep_cds_sequence_table.txt"
    paralogData = pd.read_csv(ensembl102_paralog_seq_file, header=0, sep="\t")
    paralogData = paralogData[["ensembl_gene_id","hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')
    paralogData = paralogData.iloc[:100]
    input1 = []
    for index, row in paralogData.iterrows():
        gene1 = row['ensembl_gene_id'];
        gene2 = row['hsapiens_paralog_ensembl_gene'];
        if gene1 in ensembl_uniprot_dict and gene2 in ensembl_uniprot_dict:
            gene1pdb = ensembl_uniprot_dict[gene1];
            gene2pdb = ensembl_uniprot_dict[gene2];
            input1.append((gene1,gene2,gene1pdb,gene2pdb))

    #input1= input1[:5]
    with multiprocessing.Pool(processes=thread) as pool:
        processes = pool.starmap(tmscore, input1)

    resOut = [["gene1", "gene2", "rmsd", "tm_score", "maxsub"]]
    for res in processes:
        resOut.append(res)
    df = pd.DataFrame(resOut[1:], columns=resOut[0])
    df.to_csv('Paralog_structure_similarity.txt', index=False, sep='\t')
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'TMscore jobs done!')
