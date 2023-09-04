library(ggplot2)
library(ggsignif)
##read in Key features data
Key_features <- read.csv("Key_features.tbl", sep = "\t", header = TRUE)

##calculate length of transcript
Key_features$length <- Key_features$end - Key_features$start
##make full length transcripts only
Key_features <- subset(Key_features[Key_features$feature != "exon",])
##put in ncNRA for blanks
Key_features$Type[!(Key_features$Type %in% c("coding", "ncRNA"))] <- "ncRNA"
## make expected ncRNA as own class
Key_features$class_code[is.na(Key_features$class_code)] <- "expected"
##make strucutral only ncRNA structural classcode
Key_features$class_code <- ifelse(grepl("transcript_id", Key_features$class_code), "structural", Key_features$class_code)

##make expected labelled as ncRNA
Key_features$Type[is.na(Key_features$Type)] <- "ncRNA"
## remove class codes not of interest
Key_features <- subset(Key_features, class_code %in% c('i','x', 'u', 'expected', 'tRNA', 'p', '=', "structural", "expected" ))

##plot of full-length transcript lengths

length_boxplot <- ggplot(Key_features, aes(x=Type, y=length, fill=Type)) +
  geom_boxplot() +
  scale_y_continuous(labels = function(x) x/1000, name = "Feature length (Kb)",
                     minor_breaks = seq(0, 30000, by = 1000), limits=c(0, 20000),
                     breaks = seq(0,30000, by = 5000)) +
  theme(panel.grid.major = element_line(size = 1),
        legend.position = "none") +
  geom_signif(comparisons = list(c("coding", "ncRNA")), 
              map_signif_level=TRUE)

pdf("Results/length_boxplot.pdf", width = 6, height = 6)
print(length_boxplot)
dev.off()

lengths_coding <- Key_features$length[Key_features$Type == "coding" & grepl("^[p=u]$", Key_features$class_code, ignore.case = TRUE)]
lengths_noncoding <- Key_features$length[Key_features$Type == "ncRNA"]
length_wilcox <- wilcox.test(lengths_coding, lengths_noncoding)
print(length_wilcox)


##make list of coding and non-coding genes for the python exon count script
coding_transcripts <- Key_features$Name[Key_features$Type == "coding" & grepl("^[p=u]$", Key_features$class_code, ignore.case = TRUE)]
coding_transcripts <- paste('transcript_id ', '"',coding_transcripts,'"', sep = '')
noncoding_transcripts <- Key_features$Name[Key_features$Type == "ncRNA"]
noncoding_transcripts <- paste('transcript_id ', '"', noncoding_transcripts,'"', sep = '')
write.table(coding_transcripts, "Results/coding_transcripts.txt", col.names = FALSE, row.names = FALSE, quote = FALSE )
write.table(noncoding_transcripts, "Results/noncoding_transcripts.txt", col.names = FALSE, row.names = FALSE, quote = FALSE )

##load in exon count results
exon_count_c <- read.table("Results/exon_counts_coding")
exon_count_nc <- read.table("Results/exon_counts_noncoding")

## create single dataframe with counts
exon_count_c$Type <- "mRNA"
exon_count_nc$Type <- "ncRNA"
exon_count_df <- rbind(exon_count_c,exon_count_nc)
colnames(exon_count_df) <- c("Transcript", "Count", "Type")

##plot exon counts as boxplot
exon_count_plot <- ggplot(exon_count_df, aes(x=Type, y=Count, fill=Type)) +
  geom_boxplot() +
  scale_y_continuous(name = "Exons",
                     minor_breaks = seq(0, 30000, by = 1), limits=c(0,70),
                     breaks = seq(0,70, by = 10)) +
  theme(panel.grid.major = element_line(size = 1),
        legend.position = "none")
pdf("Results/exon_count.pdf", width = 6, height = 6)
print(exon_count_plot)
dev.off()
##plot exon counts as histogram
 ggplot(exon_count_df, aes(x=Count, fill = Type)) +
  geom_histogram(binwidth = 1)



##calculate GC percentages from GC_calculate script
#GC content
GC_df_nc <- read.table("Results/Copci2021_final_noncoding.gtf_GC_count", sep = "\t")
colnames(GC_df_nc) <- c("Name", "Seq", "Percentage")
GC_df_nc$Type <- "ncRNA" 
GC_df_c <- read.table("Results/Copci2021_final_coding.gtf_GC_count", sep = "\t")
colnames(GC_df_c) <- c("Name", "Seq", "Percentage")
GC_df_c$Type <- "Coding"

GC_df <- rbind(GC_df_nc, GC_df_c)

GC_percetages_plot <- ggplot(GC_df, aes(x=Type, y=Percentage, fill = Type)) +
  geom_boxplot() +
  scale_y_continuous(name = "GC Percentage",
                     minor_breaks = seq(0, 70, by = 1), limits=c(0,70),
                     breaks = seq(0,70, by = 10)) +
  theme(panel.grid.major = element_line(size = 1),
        legend.position = "none")
pdf("Results/GC_plot.pdf", width = 6, height = 6)
print(GC_percetages_plot)
dev.off()

##ncRNA per scaffold
ggplot(per_scaffold, aes(x = Var1, y = Freq, fill = type)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(x = "Scaffold", y = "Count", fill = "transcript type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
ggsave("Results/per_scaffold.pdf")
