#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script creates a plot of the allele
freq for all the samples in the batch
"""


import vcf
import os
import csv
from do_logging import configure_logger
from collections import defaultdict
import pprint
import parse_vcf
import plot_dict
import plot_allele_freq_vcf
import itertools
import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import scipy.stats as stats
import seaborn as sns
import multiprocessing
import subprocess

logger = configure_logger("map_af")
temp_file = 'temp_af_file.csv'

def main(vcf_list):
    combine_sample_freq_dict = create_af_dict(vcf_list, logger)
    create_csv = create_giant_csv(combine_sample_freq_dict)
    create_box_plot(temp_file)
    create_violin_plot(temp_file)

def create_af_dict(vcf_list, logger):
    out_dict = defaultdict(list)
    for raw_line in vcf_list:
        vcf_file = raw_line.strip()
        sample_name = os.path.basename(vcf_file).split('_')[0]
        total_dict, snp_dict, indel_dict = parse_vcf.main(vcf_file, sample_name,logger)
        for i in total_dict.values():
            out_dict[sample_name].append(i[0][-1])
    return out_dict

def create_giant_csv(out_dict):
    with open(temp_file, 'w') as fout:
        headers = out_dict.keys()
        fout.write('\t'.join(str(i) for i in headers) + '\n')
        data = out_dict.values()
        for row in zip(*data):
            fout.write("\t".join(str(i) for i in row) + '\n')


def create_box_plot(temp_file):
    df = pd.read_csv(temp_file, sep='\t')
    fig = plt.figure(figsize=(80,25))
    ax1 = fig.add_subplot(1,1,1)
    df.boxplot()
    ax1.set(xlabel='sample', ylabel='frequency')
    ax1.yaxis.label.set_size(35)
    ax1.xaxis.label.set_size(35)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
#    plt.legend(loc='upper left')
    plt.savefig('tmp.png')
    plt.close()
    

def create_violin_plot(temp_file):
    df = pd.read_csv(temp_file, sep='\t')
    fig = plt.figure(figsize=(80,25))
    sns.set(font_scale=1.8)
    ax = sns.violinplot(df)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=14)
    ax.set(xlabel='sample', ylabel='frequency')
    plt.savefig('violin_test.png')
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='vcf_list', type=argparse.FileType('r'),
                        help='create a file with the paths to all the filtered vcf', required=True)
    args = parser.parse_args()
    main(args.vcf_list)
