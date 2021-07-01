library(data.table)
library(dplyr)
library(VennDiagram)
library(stringr)
library(ggplot2)
library(viridis)
library(gridExtra)

annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")

##split blastx
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `transcript_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `transcript_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE, keep=c(1)), by=`transcript_id`]
setnames(split.first.blastx, old=c("V1"), new=c("BlastX"))
##split blastp
blastp.results <- annotation.report[!is.na(sprot_Top_BLASTP_hit),.(sprot_Top_BLASTP_hit, `transcript_id`)]
first.blastp.hit <- blastp.results[,tstrsplit(sprot_Top_BLASTP_hit, "`", fixed = TRUE, keep=1), by = `transcript_id`]
split.first.blastp <- first.blastp.hit[,tstrsplit(V1, "^", fixed=TRUE, keep=c(1)), by=`transcript_id`]
setnames(split.first.blastp, old=c("V1"), new=c("BlastP"))
##all split annots
blast_annots_table <- full_join(split.first.blastx, split.first.blastp)
##table of genes with blastX AND blastP
blast_xp <- blast_annots_table[!is.na(blast_annots_table$BlastX) & !is.na(blast_annots_table$BlastP),]
blast_xp$gene_id <- tstrsplit(blast_xp$transcript_id, "_i", keep=c(1))
##table of transcripts with same annot from BlastX and BlastP - does still have multiple transcripts per gene though?
same_xp <- subset(blast_xp, blast_xp$BlastX==blast_xp$BlastP)

((length(same_xp$transcript_id))/(length(blast_xp$transcript_id)))*100
