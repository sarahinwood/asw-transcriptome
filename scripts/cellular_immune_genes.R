library(data.table)
library(dplyr)

##may need to change to full report once tested
trinotate_report <- fread("output/trinotate/sorted/best_annot_per_gene.csv")

##list of genes from previous papers
gene_list <- fread("data/cellular_immune_genes.csv")
##need to search trinotate file for strings in gene_list$gene_full_name
##but probably need manual checks too - some genes won't match perfectly and others may be missed
##should I be blasting against these genes somehow also?

##filter out deseq results for cellular immunity genes
ci_annots_blastx <- dplyr::filter(trinotate_report,grepl('melanization|encapsulation of foreign target|defense response to insect|melanin', gene_ontology_BLASTX))
ci_annots_blastp <- dplyr::filter(trinotate_report,grepl('melanization|encapsulation of foreign target|defense response to insect|melanin', gene_ontology_BLASTP))
ci_annots_pfam <- dplyr::filter(trinotate_report,grepl('melanization|encapsulation of foreign target|defense response to insect|melanin', gene_ontology_Pfam))

cellular_immune_GO <- full_join(ci_annots_blastx, ci_annots_blastp)

