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
output_dir="$2"
genome="$3"

# Check if the input file exists
#if [ ! -f "$input_file" ]; then
#    echo "Error: Input file not found."
#    exit 1
#fi

vcf_suff=".scored.filtered.typed.vcf"
fasta_suff=".roi_ref_kmer.fasta"

# Iterate over each file in the target directory
for family_id in "$target_dir"/*; do
    target_family="$target_dir/$family_id/"
    for target_sample in "$target_family"/*; do
        input_sample="$target_family/$target_sample"
        if [! -f "$input_file"]; then
            echo "Error: Input file not found.!"
            continue
        fi

        # Run the code for the current sample
        python get_roi_reference_kmer.py -v "$input_sample$vcf_suff" -g "$genome" -f "$input_sample$fasta_suff" -k "$input_sample.sample" -o "$output_dir/$target_sample" -p "$sample_id.roi_ref_kmer.pickle" -c 200

done
