#/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script will asses and  plot the difference is 
in total vcf
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import *
import re
import pprint
from scipy.stats import norm
import scipy.stats as stats
import logging
import shlex

def main(total_dict1, total_dict2, sample_name1, sample_name2, check_name, logger, total_dict3=False, 
         total_dict4=False, sample_name3=False, sample_name4=False):
    plot_venn = compare_vcf(total_dict1, total_dict2, total_dict3, total_dict4,
                            sample_name1, sample_name2, sample_name3, sample_name4,
                            check_name, logger)


def compare_vcf(first_dict, second_dict, third_dict, fouth_dict, 
                sample_name1, sample_name2, sample_name3, sample_name4,
                analysis_check, logger):
    sample1_sample2_outfile = "Total_missing_report_" + sample_name1 + "_"+ sample_name2 + ".txt"
    sample1_sample2_overlap, samplel_unique, sample2_unique = check_missing_in_dict(first_dict, second_dict, 
                                                                                    sample_name1, sample_name2, analysis_check, sample1_sample2_outfile)
    create_venn_plots(sample1_sample2_overlap, samplel_unique, sample2_unique, analysis_check, sample_name1, sample_name2, logger)
    if third_dict:
        sample1_sample3_outfile = "Total_missing_report_" + sample_name1 + "_"+ sample_name3 + ".txt"
        sample2_sample3_outfile = "Total_missing_report_" + sample_name2 + "_"+ sample_name3 + ".txt"
        sample1_sample3_overlap, samplel_unique, sample3_unique = check_missing_in_dict(first_dict, third_dict, sample_name1,
                                                                                        sample_name3, analysis_check, sample1_sample3_outfile)
        create_venn_plots(sample1_sample3_overlap, samplel_unique, sample3_unique, analysis_check, sample_name1, sample_name3, logger)
        sample2_sample3_overlap, sample2_unique, sample3_unique= check_missing_in_dict(second_dict, third_dict, analysis_check,
                                                                                       sample_name2, sample_name3, sample2_sample3_outfile)
        create_venn_plots(sample2_sample3_overlap, sample2_unique, sample3_unique, analysis_check, sample_name2, sample_name3, logger)
    if fouth_dict:
        sample2_sample4_outfile = "Total_missing_report_" + sample_name2 + "_"+ sample_name4 + ".txt"
        sample1_sample4_outfile = "Total_missing_report_" + sample_name1 + "_"+ sample_name4 + ".txt"
        sample3_sample4_outfile = "Total_missing_report_" + sample_name3 + "_"+ sample_name4 + ".txt"
        sample1_sample4_overlap, samplel_unique, sample4_unique = check_missing_in_dict(first_dict, fourth_dict, 
                                                                                        sample_name1, sample_name4,analysis_check,
                                                                                        sample1_sample4_outfile)
        create_venn_plots(sample1_sample4_overlap, sample1_unique, sample4_unique, analysis_check, sample_name1, sample_name4, logger)
        sample2_sample4_overlap, sample2_unique, sample4_unique = check_missing_in_dict(second_dict, fourth_dict, analysis_check,
                                                                                        sample_name2, sample_name4, sample2_sample4_outfile)
        create_venn_plots(sample2_sample4_overlap, sample2_unique, sample4_unique, analysis_check, sample_name2, sample_name4, logger)
        sample3_sample4_overlap, sample4_unique, sample4_unique = check_missing_in_dict(third_dict, fourth_dict, analysis_check,
                                                                                        sample_name3, sample_name4, sample3_sample4_outfile)
        create_venn_plots(sample3_sample4_overlap, sample3_unique, sample4_unique, analysis_check, sample_name3, sample_name4, logger)


def create_venn_plots(VCF_overlap, first_vcf_unique, second_vcf_unique, check_name, sample1, sample2, logger):
    if second_vcf_unique > 0 and first_vcf_unique > 0 and VCF_overlap > 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique,0.5, 0.5, 1.1, check_name, sample1, sample2)
    elif second_vcf_unique == 0 and first_vcf_unique == 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0, 0, 1.1, check_name, sample1, sample2)
    elif second_vcf_unique == 0 and first_vcf_unique > 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0.5, 0, 1.1, check_name, sample1, sample2)
    elif VCF_overlap == 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0.8, 0.8, 0, check_name, sample1, sample2)
    else:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0, 0.8, 1.1,check_name, sample1, sample2)


def check_missing_in_dict(some_dict1, some_dict2, sample_name1, sample_name2,analysis_check,  outfile):
    with open(analysis_check + '_variants_log.txt', 'wa+') as fout, open(outfile , 'w') as tout:
        overlap = 0
        first_vcf_unique = 0
        second_vcf_unique = 0
        first_variants = len(some_dict1)
        second_variants = len(some_dict2)
        for k, v in some_dict1.items():
            try:
                if some_dict2[k]:
                    overlap += 1
                else:
                    first_vcf_unique += 1
                    fout.write(k + '\n')
            except KeyError:
                first_vcf_unique += 1
                fout.write('Genomic' + '\t' +  str(k) + '\t' + "\t".join(str(i) for i in v).strip('[]') + '\n')
        print "overlap: " + str(overlap)
        tout.write("Total variants that overlap with {0} and {1} is {2}".format(sample_name1, sample_name2, str(overlap)) + "\n")
        print "Total variant in " + sample_name1 +  " vcf input: " + str(first_variants)
        tout.write("unique only to " + sample_name1 + " vcf input: " + str(first_vcf_unique) + '\n')
        print "unique only to " + sample_name1 +" vcf input: " + str(first_vcf_unique)
        tout.write("Total variants in "+ sample_name2 + " input: " + str(second_variants) + '\n')
        print "Total variant in " + sample_name2 + " input: " + str(second_variants)
        tout.write("unique only to " + sample_name2 + " input: " + str(second_vcf_unique) + '\n')
        print "unique only to " + sample_name2 + " " + str(second_vcf_unique)
        return overlap, first_vcf_unique, second_vcf_unique


def plot_differences(overlap, first_vcf_unique, second_vcf_unique, new_error_size, old_error_size, not_common_error_size, analysis_name, sample_name1, sample_name2):
    overlap_text = str(overlap) + " overlapping"
    new_vcf_txt = str(first_vcf_unique) + " unique " + sample_name1
    old_vcf_txt = str(second_vcf_unique) + " unique " + sample_name2
    plt.figure(figsize=(8,8))
    v = venn2(subsets = {'10': new_error_size, '01': old_error_size, '11': not_common_error_size}, set_labels = (sample_name1, sample_name2))
    v.get_patch_by_id('10').set_alpha(0.75)
    v.get_patch_by_id('10').set_color('red')
    v.get_patch_by_id('01').set_alpha(0.75)
    v.get_patch_by_id('01').set_color('orange')
    if not_common_error_size != 0:
        v.get_patch_by_id('11').set_alpha(0.75)
        v.get_patch_by_id('11').set_color('green')
        v.get_label_by_id('10').set_text(new_vcf_txt)
        v.get_label_by_id('10').set_size(10)
        v.get_label_by_id('01').set_text(old_vcf_txt)
        v.get_label_by_id('01').set_size(10)
        v.get_label_by_id('11').set_text(overlap_text)
        v.get_label_by_id('A').set_text('')
        v.get_label_by_id('B').set_text('')
        v.get_label_by_id('A').set_size(60)
        v.get_label_by_id('B').set_size(30)
        plt.annotate('1st VCF INPUT', xy = v.get_label_by_id('10').get_position(), xytext = (-30,-70), size = 'xx-large',
                     ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad=0.5', fc = 'lime', alpha = 0.3),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
        plt.annotate('2nd VCF INPUT', xy = v.get_label_by_id('01').get_position(), xytext = (30,-70), size = 'xx-large',
                     ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad = 0.5', fc = 'lime', alpha = 0.3),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad = -0.5',color = 'gray'))
        plt.savefig(analysis_name + "_VCF_VARIANT_CHECK.png")
    else:
        v.get_label_by_id('10').set_text(new_vcf_txt)
        v.get_label_by_id('10').set_size(10)
        v.get_label_by_id('01').set_text(old_vcf_txt)
        v.get_label_by_id('01').set_size(10)
        v.get_label_by_id('A').set_text('')
        v.get_label_by_id('B').set_text('')
        v.get_label_by_id('A').set_size(60)
        v.get_label_by_id('B').set_size(30)
        plt.annotate('1st VCF INPUT', xy = v.get_label_by_id('10').get_position(), xytext = (-30,-70), size = 'xx-large',
                     ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad=0.5', fc = 'lime', alpha = 0.3),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
        plt.annotate('2nd VCF INPUT', xy = v.get_label_by_id('01').get_position(), xytext = (30,-70), size = 'xx-large',
                     ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad = 0.5', fc = 'lime', alpha = 0.3),
                     arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad = -0.5',color = 'gray'))
        plt.savefig(analysis_name + "_" + sample_name1 + "_" + sample_name2 + "_VCF_VARIANT_CHECK.png")


if __name__ == "__main__":
    main()
