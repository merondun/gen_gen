#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

fasta=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/chr_W/GCA_017976375.1_bCucCan1.pri_chr_W.fa
annotation=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES/Gene_Coordinates.bed
samples=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES/BSP_Individuals.list
vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites/chr_W.AllSitesHAPMASK.vcf.gz

rm *fa
#loop through genes
for gene in $(awk '{print $4}' ${annotation}); do

	echo ${gene}
	coords=$(grep -w ${gene} ${annotation} | awk '{print $1":"$2"-"$3}')

        #loop through samples
        for sample in $(cat ${samples}); do
        echo ${sample}
        samtools faidx ${fasta} ${coords} | bcftools consensus --haplotype A ${vcf} --sample ${sample} | sed "1 s/^>.*/>${sample}/" >> ${gene}.fa

        done
done
