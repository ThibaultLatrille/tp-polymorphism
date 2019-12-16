#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import argparse
import os


class Cds(object):
    def __init__(self, chromosome, strand, name):
        """
        :chromosome : (String) Chromosome number (can also be X/Y or Z/W).
        :strand : (String) Strand on which the CDS is encoded.
        :name : (String) Name of the CDS.
        :exons : (List of 2-tuple) List of exons. Each exon is defined by a tuple (start, end),
                 where 'start' and 'end' are in absolute position in the chromosome.
        """
        self.chromosome = chromosome
        assert (strand == "+" or strand == "-")
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        if int(start_exon) <= int(end_exon):
            if self.strand == "+":
                self.exons.append((int(start_exon), int(end_exon)))
            else:
                self.exons.insert(0, (int(start_exon), int(end_exon)))
            for i in range(len(self.exons) - 1):
                if not self.exons[i][1] < self.exons[i + 1][0]:
                    print("At least one exon is overlapping with an other")

    def nt_position(self, position):
        """
        :param position: (Integer) Nucleotide position in the chromosome.
        :return: (Integer) Nucleotide position relative to the CDS.
        """
        xxxxxxxxxxxxxxxxxxxxxxxxxxxx
        return 0

    def snp_type(self, fasta_seq, position, ref_nuc, alt_nuc):
        """
        Return the SNP type given the position, reference and alternative nucleotide.
        This method also requires the fasta sequence as input.
        :param fasta_seq: (String) The fasta sequence.
        :param position: (Integer) Absolute position of the SNP in the chromosome.
        :param ref_nuc: (Character) The reference nucleotide.
        :param alt_nuc: (Character) The absolute nucleotide.
        :return: (String) The classification of the SNP, can be either "NotInCds",
                          "RefDiff", "NotIdentified", "RefStop", "Stop", "Syn" or "NonSyn".
        """
        xxxxxxxxxxxxxxxxxxxxxxxxxxx
        return "Stop"

    def empty_exon(self):
        return sum([1 for exon_len in self.exons_length() if exon_len == 1]) > 0

    def exons_length(self):
        return [j - i + 1 for i, j in self.exons]

    def seq_length(self):
        return sum(self.exons_length())


def build_dict_cds(data_path, file_name):
    print('Loading GTF file...')
    gtf_file = open("{0}/{1}".format(data_path, file_name), 'r')
    dico_cds = dict()
    for gtf_line in gtf_file:
        if gtf_line.startswith('#'):
            continue

        seqname, source, feature, start, end, score, strand, frame, comments = gtf_line.replace('\n', '').split('\t')
        if feature != 'CDS':
            continue

        transcript_find = comments.find('transcript_id')
        if transcript_find != -1 and comments.find('CCDS') != -1 and comments.find('ccds_id') != -1:
            tr_id = comments[transcript_find + 15:].split("\"")[0]
            if tr_id not in dico_cds:
                dico_cds[tr_id] = Cds(seqname, strand, tr_id)
            dico_cds[tr_id].add_exon(start, end)
    gtf_file.close()
    print('GTF file loaded.')
    return dico_cds


def most_common(lst):
    return max(set(lst), key=lst.count)


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
nucleotides = list(complement.keys())
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', required=True, type=str,
                        dest="f", metavar="<fasta>",
                        help="The relative name of the .fasta file")
    parser.add_argument('-g', '--gtf', required=True, type=str,
                        dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    parser.add_argument('-v', '--vcf', required=True, type=str,
                        dest="v", metavar="<vcf>",
                        help="The relative name of the .vcf file")
    args = parser.parse_args()

    path = os.getcwd()
    dict_cds = build_dict_cds(path, args.g)

    print('Loading fasta file...')
    dict_fasta = {}
    for fasta in SeqIO.parse(open("{0}/{1}".format(path, args.f), 'r'), 'fasta'):
        dict_fasta[fasta.id.split(".")[0]] = str(fasta.seq[:-3])
    print('Fasta file loaded.')

    stop_filename = '{0}/{1}.Stop.vcf'.format(path, args.v[:-4])
    syn_filename = '{0}/{1}.Syn.vcf'.format(path, args.v[:-4])
    nonsyn_filename = '{0}/{1}.NonSyn.vcf'.format(path, args.v[:-4])
    error_filename = '{0}/{1}.errors.tsv'.format(path, args.v[:-4])

    stop_file = open(stop_filename, 'w')
    syn_file = open(syn_filename, 'w')
    nonsyn_file = open(nonsyn_filename, 'w')
    error_file = open(error_filename, 'w')

    dict_cat_info = {"Syn": "{0} SNPs are synonymous variations",
                     "NonSyn": "{0} SNPs are non-synonymous variations",
                     "Stop": "{0} SNPs are stop variations",
                     "RefStop": "{0} SNPs have stop codon as reference amino-acid",
                     "RefDiff": "{0} SNPs retrieved from the fasta are not equal to the reference",
                     "NotIdentified": "{0} SNPs have non-identified reference or alternate amino-acid",
                     "NotInCds": "{0} SNPs are not inside the CDS"}

    dict_cat_nbr = {}
    for snp_cat in dict_cat_info:
        dict_cat_nbr[snp_cat] = 0

    vcf_file = open("{0}/{1}".format(path, args.v), 'r')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            stop_file.write(vcf_line)
            syn_file.write(vcf_line)
            nonsyn_file.write(vcf_line)
            continue

        chromo, pos, snp_id, ref, alt = vcf_line.split("\t", maxsplit=5)[:5]
        if (ref not in nucleotides) or (alt not in nucleotides):
            continue

        snp_types = []
        transcript_id_list = vcf_line[vcf_line.rfind('\t') + 1:-1].split(",")
        for transcript_id in transcript_id_list:
            snp_types.append(dict_cds[transcript_id].snp_type(dict_fasta[transcript_id], int(pos), ref, alt))

        assert len(snp_types) > 0
        if "NotInCds" in snp_types:
            dict_cat_nbr["NotInCds"] += 1
        elif "RefDiff" in snp_types:
            dict_cat_nbr["RefDiff"] += 1
        elif "NotIdentified" in snp_types:
            dict_cat_nbr["NotIdentified"] += 1
        elif "RefStop" in snp_types:
            dict_cat_nbr["RefStop"] += 1
        else:
            max_type = most_common(snp_types)
            if max_type == "Stop":
                stop_file.write(vcf_line)
                dict_cat_nbr["Stop"] += 1
            elif max_type == "Syn":
                syn_file.write(vcf_line)
                dict_cat_nbr["Syn"] += 1
            elif max_type == "NonSyn":
                nonsyn_file.write(vcf_line)
                dict_cat_nbr["NonSyn"] += 1
    vcf_file.close()

    nbr_snp_total = sum(dict_cat_nbr.values())
    error_file.write("{0} SNPs in total".format(nbr_snp_total))
    for cat, nbr in dict_cat_nbr.items():
        error_file.write("\n\n" + dict_cat_info[cat].format(nbr) + " ({0:.3f}%)".format(nbr * 100. / nbr_snp_total))

    stop_file.close()
    syn_file.close()
    nonsyn_file.close()
    error_file.close()

    print("{0} variants analyzed in total".format(nbr_snp_total))
    print("File containing {0} stop variants in {1}".format(dict_cat_nbr["Stop"], stop_filename))
    print("File containing {0} synonymous variants in {1}".format(dict_cat_nbr["Syn"], syn_filename))
    print("File containing {0} non-synonymous variants in {1}".format(dict_cat_nbr["NonSyn"], nonsyn_filename))
    nbr_errors = nbr_snp_total - dict_cat_nbr["Stop"] - dict_cat_nbr["Syn"] - dict_cat_nbr["NonSyn"]
    print("File containing {0} errors variants in {1}".format(nbr_errors, error_filename))
