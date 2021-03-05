library("data.table")
library("dplyr")
library("ggplot2")

##ASW results
asw_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")
setnames(asw_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(asw_nr_blastx, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
asw_min_evalues <- asw_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations
asw_virus <- dplyr::filter(asw_min_evalues, grepl('virus', annotation))
asw_virus <- dplyr::filter(asw_virus, !grepl('transposon', annotation))
##remove dendroctonus hit
asw_virus <- dplyr::filter(asw_virus, !grepl('Dendroctonus', annotation))
fwrite(asw_virus, "output/recip_blast/nr_blastx/viral_annots.csv")
##write list of genes vs annots
transcript_id_annot <- asw_virus[,c(1,13)]
fwrite(transcript_id_annot, "output/recip_blast/nr_blastx/viral_transcripts_annots.csv")
##use uniprot taxonomy to get genus info for each
viral_genera <- fread("output/recip_blast/nr_blastx/genus_viral_transcripts_annots.csv")
genes.per.viral.taxa <- viral_genera[,length(unique(transcript_id)), by=Genus]
fwrite(genes.per.viral.taxa ,"output/recip_blast/nr_blastx/genes_per_viral_genera.csv")
