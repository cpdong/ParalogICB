#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 10:22:59 2023

@author: cpdong
"""
import argparse, os, re;
import gzip, time;
import subprocess,shutil;
from bs4 import BeautifulSoup
from Bio import SeqIO;
import pandas as pd;
import urllib.request;
import multiprocessing;

# ~ 24*6 hrs
outdir="/gpfs/gibbs/pi/chen_sidi/cd973/0myPrj/paralogs/data/07_coevolute"
os.chdir(outdir)
thread=48;
blastp="/home/cd973/apps/ncbi-blast-2.13.0+/bin/blastp"
makeblastdb="/home/cd973/apps/ncbi-blast-2.13.0+/bin/makeblastdb"


def codeml(g1name, g1_ntseq, g1_aaseq, g2name, g2_ntseq, g2_aaseq, tmpdir): # pair-wise ks/kn calc
    #conda install -c bioconda paml
    #https://zhuanlan.zhihu.com/p/34463038

    clustalw2="/home/cd973/apps/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
    pal2nal="/home/cd973/apps/pal2nal.v14/pal2nal.pl"
    codeml="/home/cd973/miniconda3/bin/codeml"
    
    if tmpdir[-1] == '/':
        codeml_tmpdir = tmpdir + g1name + '-' + g2name
    else:
        codeml_tmpdir = tmpdir + '/' + g1name + '-' + g2name
        
    if not os.path.exists(codeml_tmpdir):
        os.makedirs(codeml_tmpdir)
    os.chdir(codeml_tmpdir)
    with open("pep.fa", "w") as f1, open("nuc.fa", "w") as f2, open("paml.tree", "w") as t1:
        f1.write(">%s\n%s\n>%s\n%s\n" % (g1name, g1_aaseq, g2name, g2_aaseq))
        f2.write(">%s\n%s\n>%s\n%s\n" % (g1name, g1_ntseq, g2name, g2_ntseq))
        t1.write("(%s, %s)\n" % (g1name, g2name))

    dNdS = 'NA'
    # clustalw2 + pal2nal
    command1 = clustalw2 + ' -ALIGN -INFILE=pep.fa -TYPE=PROTEIN -OUTFILE=clust_out.aln'
    p1 = subprocess.Popen(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    command2 = pal2nal + ' clust_out.aln nuc.fa -output paml  -nogap > paml.codon'
    p2 = subprocess.Popen(command2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];

    # codeml process
    codeml_cnt_str = "      seqfile = paml_codon_file * sequence data filename\n     treefile = paml_tree_file  * tree structure file name\n      outfile = codeml_out_file  * main result file name\n        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen\n      verbose = 1  * 0: concise; 1: detailed, 2: too much\n      runmode =-2  * 0: user tree;  1: semi-automatic;  2: automatic\n                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs\n    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n        model = 2  * models for codons:\n                   * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n                   * models for AAs or codon-translated AAs:\n                   * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n                   * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n                   * 13:3normal>0\n        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n        Mgene = 0  * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n                   * AA: 0:rates, 1:separate\n    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n        kappa = 2  * initial or fixed kappa\n    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate\n        omega = 1  * initial or fixed omega, for codons or codon-based AAs\n    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n       Malpha = 0  * different alphas for genes\n        ncatG = 4  * # of categories in dG of NSsites models\n        clock = 0  * 0: no clock, unrooted tree, 1: clock, rooted tree\n        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n       method = 0  * 0: simultaneous; 1: one branch at a time\n    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?"
    codeml_cnt_str = codeml_cnt_str.replace('paml_codon_file', codeml_tmpdir + '/paml.codon')
    codeml_cnt_str = codeml_cnt_str.replace('paml_tree_file', codeml_tmpdir + '/paml.tree')
    codeml_cnt_str = codeml_cnt_str.replace('codeml_out_file', codeml_tmpdir + '/codeml_out.txt')
    with open("codeml.cnt", "w") as t2:
        t2.write("%s\n" % codeml_cnt_str)
    p3 = subprocess.Popen(codeml + ' codeml.cnt', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    with open(codeml_tmpdir + '/codeml_out.txt', 'r') as codeml_out:
        codeml_data = codeml_out.read().splitlines();
        for line in codeml_data:
            if 'dN/dS' in line and 'dN ' in line and 'dS ' in line and '=' in line:
                dNdS = float(re.findall(r'-?\d+\.?\d*', line.split("dN/dS")[1])[0]) #re.search('dN/dS  ([\d\.\d]+)', line).group(1)

    shutil.rmtree(codeml_tmpdir)
    return [g1name, g2name, dNdS]

def blastp_self(genename, sequence, tmpdir):
    #blastp = "/Users/cpdong/Documents/Biosoft/ncbi-blast-2.13.0/bin/blastp"
    if tmpdir[-1] == '/':
        gfasta = tmpdir + genename + ".tmp.fasta"
    else:
        gfasta = tmpdir + '/' + genename + ".tmp.fasta"
    with open(gfasta, "w") as f:
        f.write(">%s\n%s\n" % (genename, sequence))
    
    out = subprocess.Popen([blastp,'-query',gfasta,'-subject',gfasta,'-outfmt','7','-evalue','1e-2'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    text = out.decode('utf-8');
    out_data = text.splitlines();
    bitscore = 'NA';
    for line in out_data:
            if not '#' in line:
                line_data = [x.strip() for x in line.split('\t') if x != ''];
                bitscore = line_data[11];
    os.remove(genename + ".tmp.fasta")
    if genelist.index(genename) % 100 ==0:
        print(genelist.index(genename))
    return [genename, float(bitscore)];

def blastp_db(genename, sequence, db_list, tmpdir):
    #blastp = "/Users/cpdong/Documents/Biosoft/ncbi-blast-2.13.0/bin/blastp"
    if tmpdir[-1] == '/':
        gfasta = tmpdir + genename + ".tmp.fasta"
    else:
        gfasta = tmpdir + '/' + genename + ".tmp.fasta"
    
    with open(gfasta, "w") as f:
        f.write(">%s\n%s\n" % (genename, sequence))
    
    if isinstance(db_list, str):
        db_list = [db_list]
    elif isinstance(db_list, list):
        db_list = db_list
    
    bitscores=[genename]; homologCheck=[genename];
    for db in db_list:
        bitscore = 0; hasHomolog = 0;
        out = subprocess.Popen([blastp,'-query',gfasta,'-db',db,'-outfmt','7','-evalue','1e-2'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];

        text = out.decode('utf-8');
        out_data = text.splitlines();
        result = [];
        for line in out_data:
            if not '#' in line:
                line_data = [x.strip() for x in line.split('\t') if x != ''];
                result.append([line_data[2],line_data[10],line_data[11]]);

        if len(result)>0:
            bitscore = float(result[0][2]) # return only the best hit
            if float(result[0][0])>35 and float(result[0][1])< 1e-04:
                hasHomolog = 1 # check whether homolog exist
        
        bitscores.append(bitscore);
        homologCheck.append(hasHomolog);
    os.remove(gfasta)

    if genelist.index(genename) % 100 ==0:
        print(genelist.index(genename))
    #return bitscores + [sum(homologCheck)]
    return [bitscores + [sum(homologCheck)], homologCheck]


if __name__ == '__main__':
    # build searching database
    ensembl102_269_taxo="https://raw.githubusercontent.com/cpdong/imParalog/refs/heads/main/data/ENSEMBL_269_species_metainfo.csv"
    taxoData = pd.read_csv(ensembl102_269_taxo, header=0, sep=",")
    taxoData = taxoData.loc[taxoData['exclude'] != 'Y']
    taxoData = taxoData[["species","taxonomy_id"]]
    taxo_dict = pd.Series(taxoData.taxonomy_id.values,index=taxoData.species).to_dict()

    # make query db from 269 species Ensembl pep fasta
    if not os.path.exists(outdir + "/db"):
        os.makedirs(outdir + "/db")

    for k,v in taxo_dict.items(): # k: species; v:taxonomy_id;
        html_page = urllib.request.urlopen('https://ftp.ensembl.org/pub/release-102/fasta/' + k + '/pep/')
        soup = BeautifulSoup(html_page)
        for link in soup.findAll('a'):
            if "pep.all.fa.gz" in link.get('href'):
                ensg_url = 'https://ftp.ensembl.org/pub/release-102/fasta/' + k + '/pep/' + link.get('href')
                if not os.path.exists(outdir + "/db/taxonomy_" + str(v)):
                    os.makedirs(outdir + "/db/taxonomy_" + str(v))
                db_prefix = outdir + "/db/taxonomy_" + str(v) + '/taxonomy_' + str(v)
                urllib.request.urlretrieve(ensg_url, "taxonomy_" + str(v) + ".fa.gz")
                with gzip.open("taxonomy_" + str(v) + ".fa.gz", 'r') as fin, open("taxonomy_" + str(v) + ".fa", "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                process = subprocess.Popen([makeblastdb,'-in',
                                            "taxonomy_" + str(v) + ".fa", '-dbtype', 'prot', 
                                            '-out', db_prefix],stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        os.remove("taxonomy_" + str(v) + ".fa.gz")
        os.remove("taxonomy_" + str(v) + ".fa")  


    ensembl102_paralog_seq_file="Ensembl102_paralog_pep_cds_sequence_table.txt"
    paralogData= pd.read_csv(ensembl102_paralog_seq_file, header=0, sep="\t")
    paralogData = paralogData.iloc[:5000]
    paras=paralogData[["geneid", "pepSeq"]].drop_duplicates(keep='first')
    mydict=pd.Series(paras.pepSeq.values,index=paras.geneid).to_dict()
    cdsData= paralogData[["geneid", "cdsSeq"]].drop_duplicates(keep='first')
    cds_dict = pd.Series(cdsData.cdsSeq.values,index=cdsData.geneid).to_dict()
    global genelist;
    genelist=paralogData.geneid.values.tolist();

    # make query db from 269 species Ensembl pep fasta
    #if not os.path.exists(outdir + "/db"):
    #    os.makedirs(outdir + "/db")
    db_list=[];
    for k,v in taxo_dict.items(): # k: species; v:taxonomy_id;
        #html_page = urllib.request.urlopen('https://ftp.ensembl.org/pub/release-102/fasta/' + k + '/pep/')
        #soup = BeautifulSoup(html_page)
        #for link in soup.findAll('a'):
            #if "pep.all.fa.gz" in link.get('href'):
            #    ensg_url = 'https://ftp.ensembl.org/pub/release-102/fasta/' + k + '/pep/' + link.get('href')
            #    if not os.path.exists(outdir + "/db/taxonomy_" + str(v)):
            #        os.makedirs(outdir + "/db/taxonomy_" + str(v))
        db_prefix = outdir + "/db/taxonomy_" + str(v) + '/taxonomy_' + str(v)
            #    urllib.request.urlretrieve(ensg_url, "taxonomy_" + str(v) + ".fa.gz")
            #    with gzip.open("taxonomy_" + str(v) + ".fa.gz", 'r') as fin, open("taxonomy_" + str(v) + ".fa", "wb") as fout:
            #        shutil.copyfileobj(fin, fout)
            #    process = subprocess.Popen([makeblastdb,'-in',
            #                                "taxonomy_" + str(v) + ".fa", '-dbtype', 'prot', 
            #                                '-out', db_prefix],stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        db_list.append(db_prefix)
                
        #os.remove("taxonomy_" + str(v) + ".fa.gz")
        #os.remove("taxonomy_" + str(v) + ".fa")

    input1 = []; input2 = []
    for k,v in mydict.items(): # k: ensgid; v:sequence;
        input1.append((k,v,db_list, outdir))
        input2.append((k,v, outdir))
        
    #input1= input1[:5]
    with multiprocessing.Pool(processes=thread) as pool:
        processes = pool.starmap(blastp_db, input1)

    resOut_1 = [['Geneid'] + [ x.split('/')[-1] for x in db_list] + ['Conservation']]
    resOut_2 = [['Geneid'] + [ x.split('/')[-1] for x in db_list]]
    for res in processes:
        resOut_1.append(res[0]);
        resOut_2.append(res[1]);
    df_homolog = pd.DataFrame(resOut_2[1:], columns=resOut_2[0]);
    df_homolog.to_csv('Paralog_homolog_conservation.txt', index=False, sep='\t'); # save a file
    df2 = pd.DataFrame(resOut_1[1:], columns=resOut_1[0]);
    
    #input2=input2[:5]
    with multiprocessing.Pool(processes=thread) as pool:
        processes = pool.starmap(blastp_self, input2)
    resOut = [['Geneid','self_bitscore']]
    for res in processes:
        resOut.append(res);
    df1 = pd.DataFrame(resOut[1:],columns=resOut[0]);
    df = pd.merge(df1, df2, on='Geneid', how='inner')
    
    urllib.request.urlretrieve("https://raw.githubusercontent.com/cpdong/public/master/ProteinHistorian_GeneAge_dollo_wagner1.0_gc36.tsv", "agefile_dollo_wagner.tsv")
    ageData=pd.read_csv("agefile_dollo_wagner.tsv", header=0, sep="\t")
    dat = pd.merge(df, ageData, on='Geneid', how='left')
    dat.to_csv('Paralog_geneAge.txt', index=False, sep='\t');
    #os.remove("agefile_dollo_wagner.tsv")
    #shutil.rmtree(outdir + '/db')
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Blast jobs done!')
    
    # Running CodeML module function from PAML
    paralog_bioMart_file="Ensembl102_pcgene_paralog_bioMart_query_result.tsv"
    para_biomart = pd.read_csv(paralog_bioMart_file, header=0, sep="\t")
    para_biomart = para_biomart[["ensembl_gene_id","hsapiens_paralog_ensembl_gene"]].drop_duplicates(keep='first')
    #para_biomart = para_biomart.head(100)
    para_biomart_list = para_biomart.values.tolist()
    input3 = []
    for item in para_biomart_list: # iterative through the paralog pairs
        gene1 = item[0]; gene2 = item[1];
        if gene1 in list(mydict.keys()) and gene2 in list(mydict.keys()):
            g1_ntseq = cds_dict[gene1]
            g1_aaseq = mydict[gene1]
            g2_ntseq = cds_dict[gene2]
            g2_aaseq = mydict[gene2]
            input3.append((gene1, g1_ntseq, g1_aaseq, gene2, g2_ntseq, g2_aaseq, outdir))

    with multiprocessing.Pool(processes=thread) as pool:
        processes = pool.starmap(codeml, input3)
    resOut = [['Gene1','Gene2','dNdS']]
    for res in processes:
        resOut.append(res);
    
    df3 = pd.DataFrame(resOut[1:],columns=resOut[0]);
    df3.to_csv('dNdS_result.txt', index=False, sep='\t');
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'CodeML jobs done!')             
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'All jobs done!')
    
