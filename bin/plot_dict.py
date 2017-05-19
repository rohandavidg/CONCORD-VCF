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

def main(total_dict1, total_dict2, check_name, logger):
    VCF_overlap, first_vcf_unique, second_vcf_unique = compare_vcf(total_dict1, total_dict2, check_name, logger)
    if second_vcf_unique > 0 and first_vcf_unique > 0 and VCF_overlap > 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique,0.5, 0.5, 1.1, check_name)
    elif second_vcf_unique == 0 and first_vcf_unique == 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0, 0, 1.1, check_name)
    elif second_vcf_unique == 0 and first_vcf_unique > 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0.5, 0, 1.1, check_name)
    elif VCF_overlap == 0:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0.8, 0.8, 0, check_name)
    else:
        get_image = plot_differences(VCF_overlap, first_vcf_unique, second_vcf_unique, 0, 0.8, 1.1,check_name)

def compare_vcf(new_dict, old_dict, analysis_check, logger):
    with open(analysis_check + '_variants_log.txt', 'w') as fout, open (analysis_check + '_First_input_variant_missing_list.txt', 'w') as vout, open(analysis_check + '_second_input_variant_missing_list.txt', 'w') as sout:
        overlap = 0
        new_vcf_unique = 0
        old_vcf_unique = 0
        new_variants = len(new_dict)
        old_variants = len(old_dict)
        for k, v in new_dict.items():
            try:
                if old_dict[k]:
                    overlap += 1
                else:
                    new_vcf_unique += 1
                    vout.write(k + '\n')
            except KeyError:
                new_vcf_unique += 1
                vout.write('Genomic' + '\t' +  str(k) + '\t' + "\t".join(str(i) for i in v).strip('[]') + '\n')

        for k, v in old_dict.items():
            try:
                if new_dict[k]:
                    pass
                else:                        
                    old_vcf_unique += 1
                    sout.write(k + '\n')
            except KeyError:
                old_vcf_unique += 1
                sout.write('Genomic' + '\t' + str(k) + '\t' + "\t".join(str(i) for i in v).strip('[]') + '\n')

        print "overlap: " + str(overlap)
        fout.write("Total variants that overlap with 1st input and 2nd input " + str(overlap) + "\n")
        print "Total variant in 1st vcf input: " + str(new_variants)
        fout.write("unique only to 1st vcf input: " + str(new_vcf_unique) + '\n')
        print "unique only to 1st vcf input: " + str(new_vcf_unique)
        fout.write("Total variants in 2nd vcf input: " + str(old_variants) + '\n')
        print "Total variant in 2nd vcf input: " + str(old_variants)
        fout.write("unique only to 2nd vcf input: " + str(old_vcf_unique) + '\n')
        print "unique only to 2nd vcf: " + str(old_vcf_unique)
        return overlap, new_vcf_unique, old_vcf_unique



def plot_differences(overlap, new_vcf_unique, old_vcf_unique, new_error_size, old_error_size, not_common_error_size, analysis_name):
    overlap_text = str(overlap) + " overlapping"
    new_vcf_txt = str(new_vcf_unique) + " unique first var"
    old_vcf_txt = str(old_vcf_unique) + " unique second var"
    plt.figure(figsize=(8,8))
    v = venn2(subsets = {'10': new_error_size, '01': old_error_size, '11': not_common_error_size}, set_labels = ('1st VCF INPUT', '2nd VCF INPUT'))
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
        plt.savefig(analysis_name + "_VCF_VARIANT_CHECK.png")


if __name__ == "__main__":
    main()
