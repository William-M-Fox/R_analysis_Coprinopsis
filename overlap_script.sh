# Iterate over each subdirectory in the current directory
for dir in */
do
    # Iterate over each GTF file in the current subdirectory
    for module_file in "$dir"/*module_*.gtf
    do
        # Create file names for ncRNA and non-ncRNA entries
        ncRNA_file="${module_file%.*}_ncRNA.gtf"
        non_ncRNA_file="${module_file%.*}_non_ncRNA.gtf"

        # Filter lines with "ncRNA" and without "ncRNA"
        grep "ncRNA" "$module_file" > "$ncRNA_file"
        grep -v "ncRNA" "$module_file" > "$non_ncRNA_file"

        # Perform bedtools intersect on the two files
        intersect_file="${module_file%.*}_overlap"
        bedtools window -w 10000 -a "$ncRNA_file" -b "$non_ncRNA_file" > "$intersect_file"
    done
done
