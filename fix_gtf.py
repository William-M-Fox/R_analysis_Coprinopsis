input_filename = "to_fix.gtf"
output_filename = "Copci2021_final_noncoding_fix.gtf"
feature_counter = {} 
with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
    for line in input_file:
        parts = line.strip().split("\t")
        if len(parts) >= 9:
            feature_type = parts[2]
            feature_info = parts[8]

            attributes = [attr.strip() for attr in feature_info.split(";")]
            transcript_id_attr = [attr for attr in attributes if attr.startswith("transcript_id")]
            
            if transcript_id_attr:
                transcript_id = transcript_id_attr[0].split('"')[1]  # Extract the transcript_id value
                
                if feature_type in feature_counter:
                    feature_counter[feature_type] += 1
                else:
                    feature_counter[feature_type] = 1

                sequential_number = feature_counter[feature_type]
                modified_transcript_id = f'{transcript_id}_{sequential_number}'

                modified_attributes = [
                    attr if not attr.startswith("transcript_id") else f'transcript_id "{modified_transcript_id}"'
                    for attr in attributes
                ]

                modified_info = "; ".join(modified_attributes)

                parts[8] = modified_info
                modified_line = "\t".join(parts)
                output_file.write(modified_line + "\n")
            else:
                output_file.write(line)  # Write unchanged lines to output
        else:
            output_file.write(line)  # Write unchanged lines to output

print("Modification complete.")