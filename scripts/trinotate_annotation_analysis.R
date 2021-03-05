library(data.table)
library(VennDiagram)
library(stringr)
library(ggplot2)
library(viridis)

#Get annotation report in right format for venn diagram
annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")
pfam <- annotation.report[!is.na(Pfam), unique(`#gene_id`)]
blastx <- annotation.report[!is.na(sprot_Top_BLASTX_hit), unique(`#gene_id`)]
kegg <- annotation.report[!is.na(Kegg), unique(`#gene_id`)]
number.genes <- annotation.report[!is.na(`#gene_id`),length(unique(`#gene_id`))]

#Draw Venn Diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")

vd <- venn.diagram(x = list("Pfam"=pfam, "BlastX"=blastx, "Kegg"=kegg), filename=NULL,
                   fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                   main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd)

#Sum of genes with any annotation
long.annotationreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam",  "eggnog", "Kegg"))
any.annotations <- long.annotationreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
any.annotations[,length(unique(`#gene_id`))]
any.annotations[,sum(any_annotations)]

#Sum of genes with any Transdecoder predicted protein
transdecoder.report <- melt(annotation.report, id.vars="#gene_id", measure.vars = c("prot_id"))
any.transdecoder <- transdecoder.report[,.(transdecoder_annotation = any(!is.na(value))),by=`#gene_id`]
any.transdecoder[,length(unique(`#gene_id`))]
any.transdecoder[,sum(transdecoder_annotation)]
sum(any.transdecoder$transdecoder_annotation==FALSE)

##transdecoder vs blastX
transdecoder <- annotation.report[!is.na(prot_id), unique(`#gene_id`)]
vd2 <- venn.diagram(x = list("Transdecoder"=transdecoder, "BlastX"=blastx), filename=NULL,
                   fill=c("#440154FF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                   main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd2)

#Annotations per genus
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
genes.per.genus <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
setkey(genes.per.genus, V1)
print(genes.per.genus)
fwrite(genes.per.genus, "output/trinotate/genes_per_genus.csv")

#meanwhile in excel sort for genus with most annotations, delete those with low no.
#I'm not interested in, and alter genus name to just genera

#plot annotations per genus - top 25 genera
plot.genes.per.genus <- fread("output/trinotate/genes_per_genus_plot.csv")
ggplot(plot.genes.per.genus, aes(x=reorder(Genus, -V1), y=V1))+
  geom_col(alpha=0.5, fill="#440154FF", colour="#440154FF")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  xlab("Genus")+ylab("Number of BlastX Annotations")

#plot annotations per class
Classes <- plot.genes.per.genus[,sum(V1), by=Class]
ggplot(Classes, aes(x=reorder(Class, -V1), y=V1))+
  geom_col(alpha=0.5, fill="#440154FF", colour="#440154FF")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  xlab("Class")+ylab("Number of BlastX Annotations")

