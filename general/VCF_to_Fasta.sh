#!/bin/bash

# Define variables
reference=$1  # Replace with the path to your reference genome in FASTA format.
vcf_file=$2  # Replace with the path to your compressed VCF file.
chromosome=$3  # Replace with the name of the chromosome of interest.
output_file=$4  # Name of the final multi-sequence FASTA file.

# Initialize output file
echo -n "" > "$output_file"

# Loop over each sample in the VCF file
for sample in $(bcftools query -l $vcf_file); do
    # Define the name of the temporary output file for each sample
    temp_output="${sample}_${chromosome}.fasta"
    
    # Generate consensus sequence for the current sample and write to a temporary output file
    bcftools consensus -f <(samtools faidx $reference $chromosome) -o "$temp_output" -s "$sample" $vcf_file
    
    # Append the temporary output file to the final multi-sequence output file
    sed "s/>.*/>${sample}/g" "$temp_output" >> "$output_file"
    
    # Remove the temporary output file
    rm "$temp_output"
done

# Print the name of the final multi-sequence FASTA file created
echo "Multi-sequence FASTA file created: $output_file"

