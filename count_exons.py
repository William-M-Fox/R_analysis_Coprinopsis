import re

import sys

if len(sys.argv) != 3:
    print("Usage: python count_exons.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]


exon_counts = {}

with open(input_file, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) >= 9 and fields[2] == "exon":
            attributes = fields[8]
            match = re.search(r'transcript_id "([^"]+)"', attributes)
            if match:
                transcript_id = match.group(1)
                exon_counts[transcript_id] = exon_counts.get(transcript_id, 0) + 1

with open(output_file, "w") as outfile:
    for transcript_id, count in exon_counts.items():
        outfile.write(f"{transcript_id}\t{count}\n")
