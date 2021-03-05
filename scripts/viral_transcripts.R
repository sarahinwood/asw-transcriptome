library(data.table)

trinotate.min.eval <- fread('output/trinotate/sorted/longest_isoform_annots.csv', na.strings = ".")
##get list of all viral genes and annots to later pull out viral transcript sequences
trinotate_virus_annots <- data.table(dplyr::filter(trinotate.min.eval, grepl('Viruses', sprot_Top_BLASTX_hit)))
fwrite(trinotate_virus_annots, "output/trinotate/viral/viral_annots.csv")
##write list of viral transcript IDs to pull out of fasta file
viral_transcripts <- trinotate_virus_annots[,2]
fwrite(list(viral_transcripts), "output/trinotate/viral/viral_transcript_ids.txt")

##virus transposon transcripts
transposon <- dplyr::filter(trinotate_virus_annots, grepl('transposon', sprot_Top_BLASTX_hit))

##sort to get viral taxa counts
trinotate_virus_annots$sprot_Top_BLASTX_hit <- tstrsplit(trinotate_virus_annots$sprot_Top_BLASTX_hit, "`", fixed=TRUE, keep=1)
trinotate_virus_annots$viral_taxa <- tstrsplit(trinotate_virus_annots$sprot_Top_BLASTX_hit, ";^", fixed=TRUE, keep=2)
genes.per.viral.taxa <- trinotate_virus_annots[,length(unique(`#gene_id`)), by=viral_taxa]
##substitute to keep last string which should be genus
#start of string, 0 or more of any character, 
genes.per.viral.taxa$viral_genera <- data.table(gsub("^.*; ", "", genes.per.viral.taxa$viral_taxa))

##have some cases where last string wasn't genus - fix these
genes.per.viral.taxa$viral_genera <- sub("Murine leukemia virus", "Gammaretrovirus", genes.per.viral.taxa$viral_genera)
genes.per.viral.taxa$viral_genera <- sub("Enterovirus G", "Enterovirus", genes.per.viral.taxa$viral_genera)
##would be good to also include whether DNA or RNA virus - search manually
fwrite(genes.per.viral.taxa, "output/trinotate/viral/genes_per_taxa_viral.csv")

