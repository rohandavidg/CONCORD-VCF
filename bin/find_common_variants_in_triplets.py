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


def main(pair1, pair2, pair3, outdir):
    analyze_trip(pair1, pair2, pair3, outdir)


def analyze_trip(pair1, pair2, pair3, outdir):
    sample1 = pair1.strip()
    sample2 = pair2.strip()
    sample3 = pair2.strip()
    sample_name1 = os.path.basename(sample1).split('_')[0]
    sample_name2 = os.path.basename(sample2).split('_')[0]
    sample_name3 = os.path.basename(sample2).split('_')[0]
    out_path = create_outdir(outdir, sample_name1, sample_name2)
    outfile = out_path + "/" + sample_name1 + '_vs_' + sample_name2 + '.report.tsv'
    sample_name1_total_dict, sample_name1_snp_dict, sample_name1_indel_dict = parse_vcf.main(sample1, 
                                                                                             sample_name1,logger)
    sample_name2_total_dict, sample_name2_snp_dict, sample_name2_indel_dict = parse_vcf.main(sample2, 
                                                                                             sample_name2,logger)
    sample_name3_total_dict, sample_name3_snp_dict, sample_name3_indel_dict = parse_vcf.main(sample3, 
                                                                                             sample_name3,logger)
    sample_name1_filtered_dict = filter_dict(sample_name1_total_dict)
    sample_name2_filtered_dict = filter_dict(sample_name2_total_dict)
    sample_name3_filtered_dict = filter_dict(sample_name3_total_dict)
    write_out_report(outfile, sample_name1_filtered_dict, sample_name1,
                     sample_name2_filtered_dict, sample_name2,
                     sample_name3_filtered_dict, sample_name3)
    common_variant_list = check_common_variant_trip(sample_name1_filtered_dict, 
                                                    sample_name2_filtered_dict, 
                                                    sample_name3_filtered_dict,
                                                    logger)
#    plot_AF_dist = plot_violin_af(common_variant_list, sample_name1, 
#                                  sample_name2, sample_name3, out_path)
#    if len(common_variant_list)> 40:
#        smaller_list = chunks(common_variant_list, 15)
#        for i, n in enumerate(smaller_list):
#            plot_AF = plot_variant_af(n, sample_name1, sample_name2, 
#                                      sample_name3, out_path, str(i))
#    else:
#        plot_AF = plot_variant_af(common_variant_list, sample_name1, 
#                                  sample_name2, sample_name3, out_path, 'freq')


def check_common_variant_trip(sample_name1_total_dict, 
                              sample_name2_total_dict, 
                              sample_name3_total_dict,
                              logger):
    common_variant = []
    set1 = set(sample_name1_total_dict)
    set2 = set(sample_name2_total_dict)
    set3 = set(sample_name3_total_dict)
    for name in set1 & set2 & set3:
        out = (name, sample_name1_total_dict[name][0][-1], 
               sample_name2_total_dict[name][0][-1], 
               sample_name3_total_dict[name][0][-1])
        common_variant.append(out)
    return common_variant


def write_out_report(outfile, sample_name1_filtered_dict, sample_name1,
                     sample_name2_filtered_dict, sample_name2,
                     sample_name3_filtered_dict, sample_name3):
    with open(outfile, 'w') as fout:
        fout.write('variant\t' + sample_name1+ '_GT\t' + sample_name1+'_ref\tref_fraction\t' + sample_name1+ '_alt\talt_fraction\t' 
                   + sample_name2 + '_GT\t' + sample_name2 + '_ref\tref_fraction\t' + sample_name2 + '_alt\talt_fraction\t' 
                   + sample_name3 + '_GT\t' + sample_name3 + '_ref\tref_fraction\t' + sample_name3 + '_alt\talt_fraction\t'
                   + 'frequency_diff\tfrequency_add\n')
        set1 = set(sample_name1_filtered_dict)
        set2 = set(sample_name2_filtered_dict)
        set3 = set(sample_name3_filtered_dict)
        for name in set1 & set2 & set3:
            #sample1
            sample1_GT = sample_name1_filtered_dict[name][0][0]
            sample1_ref = sample_name1_filtered_dict[name][0][1][0]
            sample1_ref_fraction = float(sample_name1_filtered_dict[name][0][1][0])/float(sample_name1_filtered_dict[name][0][2])
            sample1_alt = sample_name1_filtered_dict[name][0][1][1]
            sample1_alt_fraction = float(sample_name1_filtered_dict[name][0][1][1])/float(sample_name1_filtered_dict[name][0][2])
            #sample2
            sample2_GT = sample_name2_filtered_dict[name][0][0]
            sample2_ref = sample_name2_filtered_dict[name][0][1][0]
            sample2_ref_fraction = float(sample_name2_filtered_dict[name][0][1][0])/float(sample_name2_filtered_dict[name][0][2])
            sample2_alt = sample_name2_filtered_dict[name][0][1][1]
            sample2_alt_fraction = float(sample_name2_filtered_dict[name][0][1][1])/float(sample_name2_filtered_dict[name][0][2])
            #sample3
            sample3_GT = sample_name3_filtered_dict[name][0][0]
            sample3_ref = sample_name3_filtered_dict[name][0][1][0]
            sample3_ref_fraction = float(sample_name3_filtered_dict[name][0][1][0])/float(sample_name3_filtered_dict[name][0][2])
            sample3_alt = sample_name3_filtered_dict[name][0][1][1]
            sample3_alt_fraction = float(sample_name3_filtered_dict[name][0][1][1])/float(sample_name3_filtered_dict[name][0][2])
            #result
            frequency_diff = abs(sample3_alt_fraction - sample2_alt_fraction - sample1_alt_fraction)
            frequency_add = abs(sample2_alt_fraction + sample1_alt_fraction + sample3_alt_fraction)
            #outfile
            out = [name, sample1_GT, sample1_ref, sample1_ref_fraction, sample1_alt, sample1_alt_fraction, sample2_GT, sample2_ref, sample2_ref_fraction, sample2_alt, sample2_alt_fraction, sample3_GT, sample3_ref, sample3_ref_fraction, sample3_alt, sample3_alt_fraction, frequency_diff, frequency_add]
            fout.write('\t'.join(str(i) for i in out) + '\n')

if __name__ == "__main__":
    import sys
    pair1 = sys.argv[1] 
    pair2 = sys.argv[2]
    pair3 = sys.argv[3] 
    outdir = sys.argv[4]
    main(pair1, pair2, pair3, outdir)
