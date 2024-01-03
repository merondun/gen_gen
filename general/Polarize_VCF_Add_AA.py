import sys
import pysam

# Take inputs from command line arguments
outgroup_samples_file = sys.argv[1]
input_vcf_file = sys.argv[2]
output_prefix = sys.argv[3]

# Load outgroup samples from the input file
with open(outgroup_samples_file, 'r') as file:
    outgroup_samples = [line.strip() for line in file]

# Open the input VCF file
with pysam.VariantFile(input_vcf_file) as vcf_in:
    # Add 'AA' to the INFO fields of the header
    vcf_in.header.info.add('AA', '1', 'String', 'Ancestral Allele')
    # Create an output file
    with pysam.VariantFile(f"{output_prefix}.vcf.gz", 'w', header=vcf_in.header) as vcf_out:
        # Iterate over all records (SNPs)
        for record in vcf_in:
            alleles = []
            for sample in outgroup_samples:
                genotype = record.samples[sample]['GT']
                if genotype.count(None) > 1 or set(genotype) == {0, 1}:
                    # If the genotype is missing or heterozygous, mark this record as unassignable and break
                    alleles = []
                    break
                allele = record.alleles[genotype[0]]
                alleles.append(allele)
            if len(set(alleles)) > 1:
                # If more than one allele is found among the outgroup samples, assign AA=U
                ancestral_allele = 'U'
            elif len(alleles) == 0:
                # If no alleles could be determined from the outgroup samples, assign AA=U
                ancestral_allele = 'U'
            else:
                # Otherwise, the ancestral allele is the one found in the outgroup samples
                ancestral_allele = alleles[0]

            # Update the AA info
            record.info['AA'] = ancestral_allele
            vcf_out.write(record)


