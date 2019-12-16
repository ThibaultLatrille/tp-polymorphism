#!/usr/bin/env bash

RELEASE_ENSEMBL="94"

WORK_DIR=~/tp-polymorphism/data
cd ${WORK_DIR}

# Données d'annotation des séquences codantes (GTF)
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
Annotation=${WORK_DIR}/Homo_sapiens.GRCh38.98
gunzip ${Annotation}.gtf.gz

# Création du fichier interval (BED)
gtf_to_bed.py ${Annotation}.gtf
bedtools sort -i ${Annotation}.bed > ${Annotation}.sorted.bed
bedtools merge -i ${Annotation}.sorted.bed -c 4,6 -o distinct > ${Annotation}.merged.sorted.bed

### Premier Chromosome
# Données de polymorphisme (VCF)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
Polymorphisme=${WORK_DIR}/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased
gunzip ${Polymorphisme}.vcf.gz

# Filtrer le fichier VCF avec bedtools et le fichier BED
bedtools intersect -a ${Polymorphisme}.vcf -b ${Annotation}.bed -header -wb > Total.vcf

  rm ${Polymorphisme}.vcf

### Reste du Génome
for n in $(seq 2 22) "X"
do
  # Données de polymorphisme (VCF)
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${n}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
  Polymorphisme=${WORK_DIR}/ALL.chr${n}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased
  gunzip ${Polymorphisme}.vcf.gz
  
  # Filtrer le fichier VCF avec bedtools et le fichier BED
  bedtools intersect -a ${Polymorphisme}.vcf -b ${Annotation}.bed -wb >> Total.vcf
  
  rm ${Polymorphisme}.vcf
done



## Partie II, separation entre polymorphisme synonyme/non-synonyme/stop

# Données de séquences (FASTA)
wget xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
gunzip xxxxxxxxxxxxxxxxxxxxxxx

# Séparer le fichier VCF en 3 fichiers
vcf_coding_polymorphism.py xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


## Partie III, calcul de sigma2 / Va pour chacune des populations

# Recuper les informations populationnelles
wget xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Meta analyse
vcf_meta_analysis.py xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
