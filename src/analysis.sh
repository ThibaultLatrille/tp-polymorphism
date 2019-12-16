#!/usr/bin/env bash

RELEASE_ENSEMBL="94"

# Données d'annotation des séquences codantes (GTF)
wget xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
gunzip xxxxxxxxxxxxxxxxxxxxx

# Création du fichier interval (BED)
gtf_to_bed.py xxxxxxxxxxxxxxxxxxxx
bedtools sort xxxxxxxxxxxxxxxxxxxxxxxxxx
bedtools merge xxxxxxxxxxxxxxxxxxxxx

# Données de polymorphisme (VCF)
wget xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
gunzip xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Filtrer le fichier VCF avec bedtools et le fichier BED
bedtools intersect xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


## Partie II, separation entre polymorphism synonyme/non-synonyme/stop

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
