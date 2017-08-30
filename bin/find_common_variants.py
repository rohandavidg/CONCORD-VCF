#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python


"""
This script looks at common variants
in all vcf files and check the allele
frequency for cotamination
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

PYTHON='/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python'
RUN_CHECK="/dlmp/sandbox/cgslIS/rohan/cgsl_code/bin/CONCORD-VCF/bin/run_common_check.py"
QSUB = "/usr/local/biotools/oge/ge2011.11/bin/linux-x64/qsub"
LOGS_DIR="/dlmp/sandbox/cgslIS/rohan/contamination/RASFP/logs"
logger = configure_logger("compare_variants")

def main(vcf_list, outdir):
    create_sample_pairs = mutiprocess_all(vcf_list, logger, outdir)


def pair_tuple_the_list(vcf_file_list, logger, outdir):
    """
    This function create a tuple of pairs for all
    items in the list
    """
    pair_list = []
    with open(vcf_file_list) as fin:
        f = fin.readlines()
        for pair in itertools.combinations(f,2):
            pair_list.append(pair)
    return pair_list


def mutiprocess_all(vcf_file_list, logger, outdir):
    jobs = []
    pair_list = pair_tuple_the_list(vcf_file_list, logger, outdir)
    for i, n in enumerate(pair_list):
        sample_vcf1 = n[0].strip()
        sample_vcf2 = n[1].strip()
#        some_pair = [sample_vcf1, sample_vcf2]
        cmd = QSUB + " " + " -q sandbox.q " + "-N run_scan_contamination -b y" + " -l h_vmem=2G -V " + " -M gnanaolivu.rohandavid@mayo.edu" + " -e " + LOGS_DIR + " " + PYTHON + " " + RUN_CHECK + " -s1 " + sample_vcf1 + " -s2 " +  sample_vcf2 + "  -o " + outdir
        print cmd
        p =  subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate()


def analyze_pair(pair1,pair2,outdir):
    sample1 = pair1.strip()
    sample2 = pair2.strip()
    sample_name1 = os.path.basename(sample1).split('_')[0]
    sample_name2 = os.path.basename(sample2).split('_')[0]
    out_path = create_outdir(outdir, sample_name1, sample_name2)
    outfile = out_path + "/" + sample_name1 + '_vs_' + sample_name2 + '.report.tsv'
    sample_name1_total_dict, sample_name1_snp_dict, sample_name1_indel_dict = parse_vcf.main(sample1, sample_name1,logger) 
    sample_name2_total_dict, sample_name2_snp_dict, sample_name2_indel_dict = parse_vcf.main(sample2, sample_name2,logger)
    sample_name1_filtered_dict = filter_dict(sample_name1_total_dict)
    sample_name2_filtered_dict = filter_dict(sample_name2_total_dict)
    write_out_report(outfile, sample_name1_filtered_dict, sample_name1,
                     sample_name2_filtered_dict, sample_name2)
    common_variant_list =check_common_variant(sample_name1_filtered_dict, sample_name2_filtered_dict, logger)
    plot_AF_dist = plot_violin_af(common_variant_list, sample_name1, sample_name2, out_path)
    if len(common_variant_list)> 40:
        smaller_list = chunks(common_variant_list, 15)
        for i, n in enumerate(smaller_list):
            plot_AF = plot_variant_af(n, sample_name1, sample_name2, out_path, str(i))
    else:
        plot_AF = plot_variant_af(common_variant_list, sample_name1, sample_name2, out_path, 'freq')

def write_out_report(outfile, sample_name1_filtered_dict, sample_name1, 
                     sample_name2_filtered_dict, sample_name2):
    with open(outfile, 'w') as fout:
        fout.write('variant\t' + sample_name1+ '_GT\t' + sample_name1+'_ref\tref_fraction\t' + sample_name1+ '_alt\talt_fraction\t'+ sample_name2+ '_GT\t'+ sample_name2 + '_ref\tref_fraction\t'+ sample_name2 + '_alt\talt_fraction\tfrequency_diff\tfrequency_add\n')
        set1 = set(sample_name1_filtered_dict)
        set2 = set(sample_name2_filtered_dict)
        for name in set1.intersection(set2):
            sample1_GT = sample_name1_filtered_dict[name][0][0]
            sample1_ref = sample_name1_filtered_dict[name][0][1][0]
            sample1_ref_fraction = float(sample_name1_filtered_dict[name][0][1][0])/float(sample_name1_filtered_dict[name][0][2])
            sample1_alt = sample_name1_filtered_dict[name][0][1][1]
            sample1_alt_fraction = float(sample_name1_filtered_dict[name][0][1][1])/float(sample_name1_filtered_dict[name][0][2])
            sample2_GT = sample_name2_filtered_dict[name][0][0]
            sample2_ref = sample_name2_filtered_dict[name][0][1][0]
            sample2_ref_fraction = float(sample_name2_filtered_dict[name][0][1][0])/float(sample_name2_filtered_dict[name][0][2])
            sample2_alt = sample_name2_filtered_dict[name][0][1][1]
            sample2_alt_fraction = float(sample_name2_filtered_dict[name][0][1][1])/float(sample_name2_filtered_dict[name][0][2])
            frequency_diff = abs(sample2_alt_fraction - sample1_alt_fraction)
            frequency_add = abs(sample2_alt_fraction + sample1_alt_fraction)
            out = [name, sample1_GT, sample1_ref, sample1_ref_fraction, sample1_alt, sample1_alt_fraction, sample2_GT, sample2_ref, sample2_ref_fraction, sample2_alt, sample2_alt_fraction, frequency_diff, frequency_add]
            fout.write('\t'.join(str(i) for i in out) + '\n')


def filter_dict(some_dict):
    filtered_dict = {}
    for k, v in some_dict.items():
        AF = v[0][-1]
        DP = v[0][-2]
        if DP <= 10:
            pass
        else:
            filtered_dict[k] = v
    return filtered_dict


def create_outdir(outdir, sample_name1, sample_name2):
    outdir_path = outdir + '/' + sample_name1 + '_vs_' + sample_name2
    if os.path.isdir(outdir_path):
        pass
    else:
        os.makedirs(outdir_path)
    return outdir_path


def check_common_variant(sample_name1_total_dict, sample_name2_total_dict, logger):
    common_variant = []
    set1 = set(sample_name1_total_dict)
    set2 = set(sample_name2_total_dict)
    for name in set1.intersection(set2):
        out = (name, sample_name1_total_dict[name][0][-1],sample_name2_total_dict[name][0][-1])
        common_variant.append(out)
    return common_variant


def plot_variant_af(common_variant_list, sample_name1, sample_name2, outdir, iteration):
    labels = ['variant', sample_name1, sample_name2]
    df = pd.DataFrame.from_records(common_variant_list, columns=labels)
    fig = plt.figure(figsize=(30,16))
    ax1 = fig.add_subplot(1,1,1)
    x = np.linspace(0, 1, len(df['variant']))
    ax1.scatter(x, df[sample_name1], color='k', marker='x', s=550, vmin=0., vmax=0.6, label=sample_name1)
    ax1.scatter(x, df[sample_name2], color='g', marker='o', s=400, vmin=0., vmax=0.6, label=sample_name2)
    variant_label = df['variant'].tolist()
#    index = df.index.tolist()
    freq1 = df[sample_name1].tolist()
    freq2 = df[sample_name2].tolist()
    plt.xlabel("variants", fontsize=20)
    plt.ylabel("allele frequency", fontsize=20)
    for i, txt in enumerate(variant_label):
        ax1.annotate(txt, (x[i], freq1[i]), rotation=45, 
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                     arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        ax1.annotate(txt, (x[i], freq2[i]), rotation=45,
                     bbox=dict(boxstyle='round,pad=0.5', fc='green', alpha=0.5),
                     arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    else:
        pass
    plt.legend(loc='upper left')
    plt.savefig(outdir + '/' + 'AF_comparison_' + iteration + '.png')
    plt.close()


def plot_violin_af(common_variant_list, sample_name1, sample_name2, outdir):
    labels = ['variant', sample_name1, sample_name2]
    df = pd.DataFrame.from_records(common_variant_list, columns=labels)
    df = pd.DataFrame({sample_name1: df[sample_name1], sample_name2: df[sample_name2]})
    fig = plt.figure(figsize=(20,16))
    sns.set(font_scale=1.8)
    ax = sns.violinplot(df)
    ax.set(xlabel='sample', ylabel='frequency')
    plt.savefig(outdir + '/' + sample_name1 + "_" + sample_name2 +'_Allele_frequency_distribution.png')
    plt.close()

def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in xrange(0, len(l), n))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', dest='vcf_list',
                        help='create a file with the paths to all the filtered vcf', required=True)
    parser.add_argument('-o', dest='outdir',
                        help='outdir to where the output should be located', required=True)
    args = parser.parse_args()
    main(args.vcf_list, args.outdir)
