#plot number of reads mapped to a set of reference sequences
#from checkMiSeq.py

require(ggplot2)
require(RColorBrewer)
require("grid")

#palette <- c(brewer.pal(12, "Paired"), "#5C5C5C")
assign_palette <- c('phiX174'='#A6CEE3','Ecoli'='#1F78B4','Pdenitrificans'='#B2DF8A','Pacnes'='#33A02C','hg38'='#FB9A99','HCV'='#E31A1C','HIV1'='#FDBF6F','GBvirusC'='#FF7F00','GBvirusB'='#CAB2D6','HumanPegivirus2'='#6A3D9A','HBV'='#FFFF99','MMLV'='#B15928','not_referenced'='#5C5C5C')

args <- commandArgs(trailingOnly=TRUE)
print (args)

#setwd("~/Data/probedMiSeq/checkMiSeq/16_Jan_21")
#file <- '16_Jan_21_checkMiseq.txt'
#df <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
#file_prefix <- gsub(".txt", "", file)

df <- read.table(args[1], header=TRUE, sep=",", stringsAsFactors=FALSE)
file_prefix <- gsub(".txt", "", args[1])

df_4plot <- df
df_4plot$count <- as.numeric(df_4plot$count)
df_4plot$total <- as.numeric(df_4plot$total)
df_4plot$prop <- df_4plot$count/df_4plot$total * 100
df_4plot$taxon_order <- factor(df_4plot$taxon, levels = c('phiX174','Ecoli','Pdenitrificans','Pacnes','hg38','HCV','HIV1','GBvirusC','GBvirusB','HumanPegivirus2','HBV','MMLV','not_referenced'))

df_4plot$sample_num <- as.numeric(gsub("S", "", df_4plot$snum))
df_4plot$sample_labels <- paste(df_4plot$sample,df_4plot$snum, sep = '_')

df_4plot$sample_labels <- factor(df_4plot$sample_labels, levels=rev(df_4plot[order(df_4plot$sample_num), "sample_labels"]))

all_prop <- (ggplot(df_4plot, aes(x = sample_labels, y=prop, fill=taxon_order, order=taxon_order)) + 
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual(values = assign_palette) +
  ylab("proportion of total reads hitting references") +
  coord_flip() +
  ggtitle(paste(file_prefix, "\nproportion of total reads hitting references (%)")) +
  theme(plot.title = element_text(lineheight=.8, face="bold"), axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=12, color = "black"), axis.title.x = element_text(size=14), axis.text.y = element_text(color = "black", size=12), axis.title.y = element_blank(), legend.title = element_blank(), legend.key.size = unit(0.7, "cm"), legend.text = element_text(size=12)))

#all_prop
filename = paste(file_prefix, 'checkMiSeq.pdf', sep='_')
ggsave(filename, plot = all_prop)


