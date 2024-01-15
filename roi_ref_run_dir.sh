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
        input_dir="$target_family/$target_sample/"
        #if [! -f "$input_file"]; then
        #    echo "Error: Input file not found.!"
        #    continue
        #fi

        # Run the code for the current sample
        python get_roi_reference_kmer.py -v "$input_dir$target_sample$vcf_suff" -g "$genome" -f "$input_dir$target_sample$fasta_suff" -k "$input_dir$target_sample.sample" -o "$input_dir" -p "$input_dir$target_sample.roi_ref_kmer.pickle" -c 200
        kmc -v -fa -k31 -ci1 -t10 -m8 "$input_dir$target_sample$fasta_suff" "$input_dir$target_sample.roi_ref_kmer" "input_dir"
        kmc_tools simple "$input_dir$target_sample.sample" "$input_dr$target_sample.roi_ref_kmer" intersect "$input_dir$target_sample.roi_ref_kmer_count" -ocleft
        kmc_tools transform "$input_dir$target_sample.roi_ref_kmer_count" dump "$input_dir$target_sample.roi_ref_kmer_count_dump.txt"
        python save_kmerdb_pickle.py -i "$input_dir$target_sample.roi_ref_kmer_count_dump.txt" -p "$input_dir$target_sample.roi_ref_kmer_count.pickle"
        echo "kmer db pickle object saved at $input_dir$target_sample.roi_ref_kmer_count.pickle"

done
