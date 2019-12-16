#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

RED = "#EB6231"
BLUE = "#5D80B4"
GREEN = "#8FB03E"

dico_name = {"CHB": "Han Chinese in Beijing, China",
             "JPT": "Japanese in Tokyo, Japan",
             "CHS": "Southern Han Chinese",
             "CDX": "Chinese Dai in Xishuangbanna, China",
             "KHV": "Kinh in Ho Chi Minh City, Vietnam",
             "CEU": "Utah Residents (CEPH) with European Ancestry",
             "TSI": "Toscani in Italia",
             "FIN": "Finnish in Finland",
             "GBR": "British in England and Scotland",
             "IBS": "Iberian Population in Spain",
             "YRI": "Yoruba in Ibadan, Nigeria",
             "LWK": "Luhya in Webuye, Kenya",
             "GWD": "Gambian in Western Divisions in the Gambia",
             "MSL": "Mende in Sierra Leone",
             "ESN": "Esan in Nigeria",
             "ASW": "Americans of African Ancestry in SW USA",
             "ACB": "African Caribbeans in Barbados",
             "MXL": "Mexican Ancestry from Los Angeles USA",
             "PUR": "Puerto Ricans from Puerto Rico",
             "CLM": "Colombians from Medellin, Colombia",
             "PEL": "Peruvians from Lima, Peru",
             "GIH": "Gujarati Indian from Houston, Texas",
             "PJL": "Punjabi from Lahore, Pakistan",
             "BEB": "Bengali from Bangladesh",
             "STU": "Sri Lankan Tamil from the UK",
             "ITU": "Indian Telugu from the UK",
             "AFR": "African",
             "AMR": "Admixed American",
             "EAS": "East Asian",
             "EUR": "European",
             "SAS": "South Asian"}


def count_comments(filename):
    """Count comment lines (those that start with "##")
    :param filename: (String) A file.
    :return: (Integer) the .
    """
    nbr_comments = 0
    with open(filename) as fh:
        for line in fh:
            if line.startswith('##'):
                nbr_comments += 1
            else:
                break
    return nbr_comments + 1, line.split('\t') + ['chr', 's', 'e', 'tr']


def epistasis(array):
    """
    Compute the sigma2 / Va.
    :param array: (2-d array of float) SNPs are in rows, individuals in columns.
    :return: (Float) Epistasis = sigma2 / Va.
    """
    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    return 1.0


def bootstrap(nbr_variants, array, bootstrap_repetition=1000):
    """
    Extract random samples of and compute epistasis for each sample.
    :param nbr_variants: (Integer) The number of SNPs in each random sample.
    :param array: (2-d array of float) SNPs are in rows, individuals in columns.
    :param bootstrap_repetition: (Integer) The number of random sample to draw.
    :return: (List of float) Epistasis computed for each random sample.
    """
    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    return [1.0] * bootstrap_repetition


def p_value(bootstrap_values, value):
    """
    Return the p-value.
    :param bootstrap_values: (List of float) The number of SNPs in each random sample.
    :param value: (2-d array of float) SNPs are in rows, individuals in columns.
    :return: (Float) p-value between 0 and 1.
    """
    return len([1 for x in bootstrap_values if x < value]) / len(bootstrap_values)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--stop', required=True, type=str,
                        dest="stop", metavar="<vcf>",
                        help="The relative name of the stop .vcf file")
    parser.add_argument('--syn', required=True, type=str,
                        dest="syn", metavar="<vcf>",
                        help="The relative name of the synonymous .vcf file")
    parser.add_argument('--nonsyn', required=True, type=str,
                        dest="nonsyn", metavar="<vcf>",
                        help="The relative name of the non-synonymous .vcf file")
    parser.add_argument('--panel', required=True, type=str,
                        dest="panel", metavar="<panel>",
                        help="The relative name of the .panel file")
    parser.add_argument('--cutoff', required=False, default=0.01, type=float,
                        dest="cutoff", metavar="<cutoff>",
                        help="The cut-off for minor allele frequency.")

    args = parser.parse_args()

    tsv_file = open("{0}/meta_analysis_{1}.tsv".format(os.getcwd(), args.cutoff), 'w')
    tsv_file.write(
        '\t'.join(['Population', 'NbrIndividuals', 'NbrStops', 's^2/Va', 'Pvalue', 'Significant']) + '\n')

    panel = pd.read_csv("{0}/{1}".format(os.getcwd(), args.panel), sep='\t', usecols=('sample', 'pop'))

    genotypes = dict()
    for snp_type in ['stop', 'syn', 'nonsyn']:
        filepath = "{0}/{1}".format(os.getcwd(), getattr(args, snp_type))
        nb, header = count_comments(filepath)
        individuals = set(panel['sample']).intersection(header)
        genotypes[snp_type] = pd.read_csv(filepath, sep='\t', skiprows=nb, names=header,
                                          converters={k: lambda x: x.count("1") for k in individuals},
                                          usecols=individuals)

    for pop in set(panel['pop']):
        individuals = panel[panel["pop"] == pop]["sample"]
        nb_alleles = 2 * len(individuals)
        print("{0} population with {1} individuals".format(pop, len(individuals)))

        filtered_array = dict()
        for snp_type in ['stop', 'syn', 'nonsyn']:
            df = genotypes[snp_type].filter(items=individuals)
            df_sum = df.sum(axis=1)
            filt_freq = ((0 < df_sum) & (df_sum < nb_alleles * args.cutoff)) | (
                    (nb_alleles * (1 - args.cutoff) < df_sum) & (df_sum < nb_alleles))
            filtered_array[snp_type] = df[filt_freq].values

        stop_epistasis = epistasis(filtered_array['stop'])
        nbr_stops = len(filtered_array['stop'])
        if nbr_stops > 1:
            syn_bootstrap = bootstrap(nbr_stops, filtered_array['syn'])
            non_syn_bootstrap = bootstrap(nbr_stops, filtered_array['nonsyn'])

            my_dpi = 256
            fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            p_val = p_value(syn_bootstrap, stop_epistasis)

            result_line = [pop, str(len(individuals)), str(nbr_stops), str(stop_epistasis), str(p_val)]
            if p_val < 0.05:
                result_line += ["True"]
                print("\t p-value={0:3g}, the statistical test is significant ".format(p_val))
            else:
                result_line += ["False"]
                print("\t p-value={0:3g}, the statistical test is not significant ".format(p_val))
            tsv_file.write('\t'.join(result_line) + '\n')

            plt.title('{0} {1},\n p-value={2:3g} ({3} alleles up to a minor allele count of {4})'.format(
                len(individuals), dico_name[pop], p_val, nbr_stops, args.cutoff), fontsize=8)
            bins = 50
            syn_hist, _, _ = plt.hist(syn_bootstrap, bins, density=1, facecolor=BLUE,
                                      alpha=0.4, label='Synonymous (10,000 resampling)')
            non_syn_hist, _, _ = plt.hist(non_syn_bootstrap, bins, density=1, facecolor=GREEN,
                                          alpha=0.4, label='Non-Synonymous (10,000 resampling)')
            y_max = 1.2 * max((max(syn_hist), max(non_syn_hist)))
            plt.ylim((0, y_max))
            plt.plot((stop_epistasis, stop_epistasis), (0, y_max), linewidth=3, color=RED,
                     label=(r'LoF ($\sigma^2/V_{A}$' + '={0:3g})'.format(stop_epistasis)))
            plt.xlabel(r'$\sigma^2/V_{A}$', fontsize=8)
            plt.ylabel('Density', fontsize=8)
            plt.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig("{0}/analysis_{1}_{2}.svg".format(os.getcwd(), pop, args.cutoff), format="svg")
            plt.close()
        else:
            print('\tNo stop variants for this population')

    tsv_file.close()
    print("Analysis completed")
