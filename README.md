make_gtf.R
- uses a gffcomapre output as a base gtf
- combines a base gtf with InterProScan results, RNAcode, PLEK, and CPC2 coding calcualtor results to create afinal gtf of coding and non-coding sequences

Feature_analysis.R 
- Uses the Key_features output from make_gtf.R to graph transcript length, GC percentage and exon count for mRNA and ncRNA and calcualtes ncRNA per scaffold

Deseq_analysis.R
- Uses DeSeq2 to perform differential expression analysis, graph heatmaps for ncRNA and mRNA seperately.
- Graph a 100% stacked column graph of transcript type and direction for each stage, stage specificty and PCA of samples

count_exons.py
Usage: python count_exons.py input_file output_file
- takes a gtf input and calculates exons per transcript

fix_gtf.py
- add sequential numbers to structural ncRNA to differentiate them

create_gtf_files.sh
- takes network correlation moduels as input to make gtf files for each module
- 
overlap script.sh
- takes a module gtf file as input and seperates coding and non-coding transcripts 
- uses bedtools windows to identify genes within 10kb to identfy potential interaction  (ncRNA-mRNA only)

combine_files.sh
- combines overlap files for each sample
