import pandas as pd;
import numpy as np;
import os, gzip;
import subprocess,shutil;
from Bio import SeqIO;
from multiprocessing import Pool;

os.chdir("/path/to/diretory")
fasta = 'gencode.v36.pc_translations.fa.gz';
biomart_file = 'Ensembl102_pcgene_paralog_bioMart_query_result.tsv';

def extractlongestFasta(input_file):
    if not os.path.exists('tmp'):
        os.makedirs('tmp')
    txlist=[];
    if input_file.lower().endswith('.gz'):
        infile=gzip.open(input_file,"rt");
    elif input_file.lower().endswith(('.fasta', '.fa')):
        infile=open(input_file,"rt");
    else:
        infile=None;
    for record in SeqIO.parse(infile,"fasta"):
          ensg=[x for x in str(record.id).split('|') if "ENSG" in x ][0][:15]
          enst=[x for x in str(record.id).split('|') if "ENST" in x ][0]
          symb=str(record.id).split('|')[-2]
          txlist.append([ensg,symb,len(str(record.seq)),enst,str(record.seq)])

    glist=list(set([x[0] for x in txlist]))
    print("Total pcgene", len(glist))
    for g in glist:
        gtxlist=[x for x in txlist if x[0]==g]
        max_len = max([x[2] for x in gtxlist])
        gItem=[x for x in gtxlist if x[2] == max_len][0]
        f1.write('|'.join(gItem[:2]) + '|' + str(gItem[2]) + '\t' + gItem[-1] + '\n')
        with open('./tmp/' + g + '.fasta', 'w') as fout:
            fout.write('>' + '|'.join(gItem[:2]) + '|' + str(gItem[2]) + '\n' + gItem[-1])

def blastp2(gene1, gene2):
    gene1_fasta = './tmp/' + gene1 + '.fasta'
    gene2_fasta = './tmp/' + gene2 + '.fasta'
    out = subprocess.Popen(['blastp','-query', gene1_fasta, '-subject', gene2_fasta, '-outfmt', '7', '-evalue', '1e-0'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    text = out.decode('utf-8');
    out_data = text.splitlines();
    result = [];
    for line in out_data:
        if not '#' in line:
            line_data = [x.strip() for x in line.split('\t') if x != ''];
            result.append(line_data);
    if len(result)>1:
        pval = [float(x[-2]) for x in result];
        return result[pval.index(min(pval))]
    elif len(result)==1:
        return result[0]
    else:
        return [gene1, gene2] + ['']*10;

def combine_blastp(biomart_list):
    if biomart_list and len(biomart_list)>0:
        gene1 = biomart_list[0][:15]
        gene2 = biomart_list[1][:15]
        # load blastp result files
        blastp_data = blastp2(gene1, gene2);
        biomart_list = biomart_list + list(np.array(blastp_data)[[2,3,4,5,10,11]]);
        biomart_list[:2] = blastp_data[:2]
        return biomart_list;
    else:
        return None;

if __name__ == '__main__':

    extractlongestFasta(fasta, './');
    unigene = [x[:15] for x in os.listdir('./tmp')]
    print(len(unigene))

    biomart = [];
    with open(biomart_file) as f1: 
        next(f1)
        for line in f1:
            line_data = line.rstrip('\n').split("\t");
            if line_data[0] in unigene and line_data[1] in unigene:
                biomart.append(list(np.array(line_data)[[0, 1, 8, 9, 10, 11]]));
    print(len(biomart))

    with Pool(processes=thread) as pool:
        process1 = pool.map(combine_blastp, biomart); #multiple process
    
    resOut = [["query_gene","paralog_gene","hsapiens_paralog_subtype","hsapiens_paralog_orthology_type","hsapiens_paralog_perc_id",
            "hsapiens_paralog_perc_id_r1", "% identity","alignment length","mismatches","gap opens","evalue","bit.score"]]
    for res in process1:
        resOut.append(res);

    df = pd.DataFrame(resOut[1:],columns=resOut[0]);
    df.to_csv(outdir + 'Paralog_pairs_sequence_identify_BLASTP.txt', index=False, sep='\t');
    shutil.rmtree('./tmp')
