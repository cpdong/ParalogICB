library("biomaRt")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=102)
attrlist=listAttributes(ensembl)$name
para.attr<- c("ensembl_gene_id", attrlist[grepl("paralog", attrlist)])

inputlist<- getBM(attributes = "ensembl_gene_id", filters = c("transcript_biotype","with_hsapiens_paralog"),
                  values = list(transcript_biotype="protein_coding", with_hsapiens_paralog=TRUE), mart = ensembl)$ensembl_gene_id
result<- getBM(attributes = para.attr, filters = "ensembl_gene_id", values = inputlist, mart = ensembl)
write.table(result, "Ensembl102_pcgene_paralog_bioMart_query_result.tsv", row.names=F, quote=F, sep="\t")
