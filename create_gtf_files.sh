gtf_file="Copci2021_final_trans.gtf"

# Iterate over each subdirectory in the current directory
for dir in */
do
    # Iterate over each file in the current subdirectory
    for file in "$dir"/*
    do
        # Check if the file is a text file
        if [[ "$file" == *.txt ]]
        then
            # Extract transcript IDs from the file
            transcript_ids=$(grep -o -F -f "$file" "$gtf_file")

            # Create a new GTF file using the extracted transcript IDs
            output_file="${file%.txt}.gtf"
            grep -F "$transcript_ids" "$gtf_file" > "$output_file"

            echo "Created $output_file"
        fi
    done
done