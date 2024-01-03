#!/bin/bash

usage() {
    echo "Usage: $0 <GFF_FILE> <GENOME_FILE> <BED_FILE>"
    echo "  GFF_FILE: Path to the GFF file."
    echo "  GENOME_FILE: Path to the genome file containing chromosome lengths."
    echo "  BED_FILE: Path for the output BED file."
    echo "Example: $0 my_annotations.gff my_genome.txt my_output.bed"
}

if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

# Define your GFF file and the output BED file
GFF_FILE=$1
GENOME_FILE=$2
BED_FILE=$3

# Define N KB for promoter region (e.g., 1500 for 1.5KB)
PROMOTER_REGION=1500

# Extract exons, CDS, and genes
awk '$3=="exon"' $GFF_FILE > exons.gff
awk '$3=="gene"' $GFF_FILE > genes.gff

# Convert GFF to BED
gff2bed < exons.gff | bedtools sort -i - | awk '{OFS="\t"}{print $1, $2, $3, $8, ".", $6}' | uniq > exons.bed
gff2bed < genes.gff | bedtools sort -i - | awk '{OFS="\t"}{print $1, $2, $3, $8, ".", $6}' | uniq > genes.bed

# Find introns per gene
bedtools subtract -a genes.bed -b exons.bed | awk '{OFS="\t"}{print $1, $2, $3, "intron", $5, $6}' > introns.bed

# Find promoters (regions N bases upstream of CDS)
bedtools flank -i genes.bed -g $GENOME_FILE -l $PROMOTER_REGION -r 0 -s | awk '{OFS="\t"}{print $1, $2, $3, "promoter", $5, $6}' > promoters.bed
cp promoters.bed $BED_FILE

# Exclude promoter regions from exons and add remaining exons
bedtools subtract -a exons.bed -b promoters.bed | bedtools sort -i - >> $BED_FILE
bedtools sort -i $BED_FILE > tmp; mv tmp $BED_FILE

# Exclude promoter and exon regions from introns and add remaining introns
bedtools subtract -a introns.bed -b $BED_FILE | bedtools sort -i - >> $BED_FILE
bedtools sort -i $BED_FILE > tmp; mv tmp $BED_FILE

# Exclude promoter and exon and intron regions from intergenic and add remaining intergenic
awk '{OFS="\t"}{print $1, 0, $2}' $GENOME_FILE | bedtools subtract -a - -b $BED_FILE | bedtools sort -i - | awk '{OFS="\t"}{print $1, $2, $3, "intergenic", ".", "+"}' >> $BED_FILE
bedtools sort -i $BED_FILE > tmp; mv tmp $BED_FILE

awk '{OFS="\t"}{print $1, $2, $3-1, $4, $5, $6}' $BED_FILE | bedtools merge -i - -c 4 -o distinct > tmp1
awk -F'\t' '{
    split($4, a, ",");
    for (i in a) {
        if (a[i] == "promoter") { feature = "promoter"; break; }
        else if (a[i] == "exon") { feature = "exon"; }
        else if (a[i] == "intron") { feature = "intron"; }
        else if (a[i] == "intergenic") { feature = "intergenic"; }
    }
    print $1, $2, $3+1, feature
}' OFS='\t' tmp1 > $BED_FILE

#total size annotated
size=$(awk '{print $3-$2}' $BED_FILE | datamash sum 1)
raw=$(awk '{print $2}' $GENOME_FILE | datamash sum 1)
echo "RAW: ${raw} and ANNOTATED: ${size}"
awk '{print $4, $3 - $2}' $BED_FILE | awk '{OFS="\t"}{print $1, $2}' | datamash -s -g 1 sum 2

rm tmp1 introns.bed exons.bed genes.bed promoters.bed exons.gff genes.gff 
