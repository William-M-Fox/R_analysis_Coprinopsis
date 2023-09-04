for dir in */; do
    subdir="${dir%/}"  # Remove trailing slash (/) from directory name
    combined_file="${subdir}.combined.txt"

    # Combine *.overlap files into a single file
    cat "$subdir"/*_overlap > "$combined_file"
done