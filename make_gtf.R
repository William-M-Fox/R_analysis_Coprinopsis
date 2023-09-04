library(reshape2)
library(tidyverse)
setwd("/Volumes/student_users/williamfox/R_paper/R_paper")
#load in gffcompare output with trnascripts_ids as a seperate column
Copci2021_gtf <- read.csv("CC_ALL_Copci2021.annotated_known_all_structural_all_ids.gtf", sep = "\t", header = FALSE, quote = "")
colnames(Copci2021_gtf) <-  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "Name" )
Copci2021_gtf$Name <- ifelse(Copci2021_gtf$Name == "", Copci2021_gtf$attribute, Copci2021_gtf$Name)
Copci2021_gtf$Name <- gsub('"transcript_id"', '', Copci2021_gtf$Name)

# load in coding predictions
hypothetical_prediction_df <- read.delim("hypothetical_Prediction", header =FALSE)
novel_prediction_df <- read.delim("new_novel_prediction", header = FALSE)
colnames(hypothetical_prediction_df) <- c("Program", "Prediction", "Sequence")
colnames(novel_prediction_df) <- c("Program", "Prediction", "Sequence")

# Reformat prediction
hypothetical_prediction_df_wide <- dcast(hypothetical_prediction_df, Sequence ~ Program, value.var = "Prediction")
novel_prediction_df_wide = dcast(novel_prediction_df, Sequence ~ Program, value.var = "Prediction")

# replace "na" from RNAcode column
hypothetical_prediction_df_wide$RNAcode <- hypothetical_prediction_df_wide$RNAcode %>% replace_na('None')
novel_prediction_df_wide$RNAcode <- novel_prediction_df_wide$RNAcode %>% replace_na('None')

## identify novel identified as noncoding
novel_noncoding_all <- novel_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "Noncoding")
novel_noncoding_none <- novel_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "None")
novel_ncRNA_prediction <- rbind(novel_noncoding_all, novel_noncoding_none)

##identify novel identified as coding
novel_coding <- novel_prediction_df_wide %>% filter(CPC2 == "Coding" | PLEK == "Coding" | RNAcode == "Coding")

#### identify hypothetical identified as noncoding
hypothetical_noncoding_all <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "Noncoding")
hypothetical_noncoding_none <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "None")
hypothetical_ncRNA_prediction <- rbind(hypothetical_noncoding_all, hypothetical_noncoding_none)

##identify hypothetical identified as coding
hypothetical_coding <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Coding" | PLEK == "Coding" | RNAcode == "Coding")

#combining coding predictions
coding_putative <- rbind(hypothetical_coding, novel_coding)
noncoding_putative<- rbind(hypothetical_ncRNA_prediction, novel_ncRNA_prediction)
all_putative_analysis <- rbind(noncoding_putative,coding_putative)


#make into a format comparable to the gtf
noncoding_putative<- rbind(hypothetical_ncRNA_prediction, novel_ncRNA_prediction)
noncoding_putative <- separate(noncoding_putative, col = Sequence, into = c("seqname", "start", "end"), sep = "[:-]")
noncoding_putative$start <- as.numeric(noncoding_putative$start) + 1
noncoding_putative$end <- as.numeric(noncoding_putative$end)
noncoding_putative$seqname <- gsub(">","",noncoding_putative$seqname)
noncoding_putative[,7] <- "ncRNA"

#combine coding results to gtf
Copci2021_noncoding_gtf <- left_join(Copci2021_gtf, noncoding_putative, by= c("seqname", "start", "end"), multiple= "first")

# class codes into a seperate column
Copci2021_noncoding_gtf$class_code[Copci2021_noncoding_gtf$feature == "transcript"] <- gsub('.*class_code "([=a-z]).*','\\1', Copci2021_noncoding_gtf$attribute[Copci2021_noncoding_gtf$feature == "transcript"])

#load hypothetical interpro 
hypothetical_interpro <- read.csv("hypothetical_interpro.tabular", sep = "\t", header=FALSE)
colnames(hypothetical_interpro) [1] = "Name"
hypothetical_interpro <- hypothetical_interpro[grep("Pfam", hypothetical_interpro$V4),]

#add to gtf
Copci2021_noncoding_interpro <- left_join(Copci2021_noncoding_gtf, hypothetical_interpro, by="Name", multiple = "first")
Copci2021_noncoding_interpro <- Copci2021_noncoding_interpro[,c(1:10,14:15,18:20, 26:28)]



#add novel_interpro
novel_interpro <- read.csv("novel_interpro.tabular", sep = "\t", header=FALSE)
colnames(novel_interpro) [1] = "Name"
novel_interpro <- novel_interpro[grep("Pfam", novel_interpro$V5),]
novel_interpro <- separate(novel_interpro, col = Name, into = c("seqname", "start", "end"), sep = "[:-]")
novel_interpro$start <- as.numeric(novel_interpro$start)
novel_interpro$end <- as.numeric(novel_interpro$end)
Copci2021_noncoding_interpro <- left_join(Copci2021_noncoding_interpro, novel_interpro, by= c("seqname", "start", "end"), multiple = "first")
Copci2021_noncoding_interpro$V7.x <- ifelse(
  is.na(Copci2021_noncoding_interpro$V7.x) & 
    Copci2021_noncoding_interpro$feature == "transcript",
  "coding",
  Copci2021_noncoding_interpro$V7.x
)
Copci2021_noncoding_interpro$class_code[grepl("unassigned_transcript", Copci2021_noncoding_interpro$Name) & Copci2021_noncoding_interpro$feature == "transcript"] <- "tRNA"
Copci2021_noncoding_interpro$V7.x[Copci2021_noncoding_interpro$class_code == "tRNA"] <- "ncRNA"

colnames(Copci2021_noncoding_interpro)[colnames(Copci2021_noncoding_interpro) == "V7.x"] <- "Type"

Copci2021_noncoding_interpro$Type[Copci2021_noncoding_interpro$V4.x == "Pfam" | Copci2021_noncoding_interpro$V5.y == "Pfam"] <- "coding"


#Make a dataframe of key features for analysis
Key_features <- Copci2021_noncoding_interpro[,c(1:5, 10:12)]
write.table(Key_features, "Key_features.tbl", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE )

#table of Go terms for each transcript
go_terms <- Copci2021_noncoding_interpro[, c("Name", "V14.x", "V15")]
go_terms$V14.x <- replace(go_terms$V14.x, go_terms$V14.x == "-", NA)
go_terms$V15 <- replace(go_terms$V15, go_terms$V15 == "-", NA)
go_terms <- go_terms %>%
  filter(!is.na(V14.x) | !is.na(V15))
go_terms <- go_terms %>%
  mutate(V14.x = coalesce(V14.x, V15),
         V15 = ifelse(is.na(V14.x), V15, NA))
go_terms <- go_terms[, c("Name", "V14.x")]


#clean up table
Copci2021_noncoding_interpro <- Copci2021_noncoding_interpro[,c(1:18,22:24,30:32)]

# make table into gtf style 
Copci2021_noncoding_interpro$combined_attributes <- apply(Copci2021_noncoding_interpro[, c(9, 11:24)], 1, function(x) {
  paste(ifelse(is.na(x) | x == "", "", x), collapse="; ")
})

Copci2021_noncoding_interpro  <- subset(Copci2021_noncoding_interpro, select=c(1:8, 25))
Copci2021_noncoding_interpro$combined_attributes <- gsub("; ;", "", Copci2021_noncoding_interpro$combined_attributes)
write.table(Copci2021_noncoding_interpro, "Results/Copci2021_final.gtf", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
