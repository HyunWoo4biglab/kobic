#!/bin/bash

# Input file containing sample information
#input_file="sample_list.txt"

# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Input file containing sample information
input_file="$1"
outuput_dir="$2"
genome="$3"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file not found."
    exit 1
fi

vcf_suff=".scored.filtered.typed.vcf"
fasta_suff=".roi_ref.fasta"


# Iterate over each line in the input file
while IFS= read -r line; do
    # Extract sample id and input file path from the line
    sample_id=$(echo "$line" | awk '{print $1}')
    input_path=$(echo "$line" | awk '{print $2}')

    # Check if the input file for the sample exists
    if [ ! -f "$input_path" ]; then
        echo "Error: Input file not found for sample $sample_id."
        continue
    fi

    # Run the python code for the current sample
    python get_roi_reference_kmer.py -i "$input_path"
    python get_roi_reference_kmer.py -v "$sample_id$vcf_suff" -g "$genome" -f "$sample_id$fasta_suff" -k "$sample_id.sample" -o "$outuput_dir/$sample_id" -p "$sample_id.roi_ref_kmer.pickle" -c 200


done < "$input_file"
