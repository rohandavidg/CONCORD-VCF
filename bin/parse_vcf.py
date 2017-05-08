#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

""" 
parse a vcf file 
tailored for CLC vcf for the moment
"""

import vcf
import os
import csv
from logging import logging
from collections import defaultdict
import pprint

def main(vcf_file):
    vcf_snp_dict, vcf_indel_dict, vcf_del_dict, vcf_transition_dict = vcf_to_dict(vcf_file)
    return vcf_snp_dict, vcf_indel_dict, vcf_del_dict, vcf_transition_dict


def vcf_to_dict(vcf_file):
    snp_dict = {}
    indel_dict = {}
    del_dict = {}
    transition_dict = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.is_snp:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'], 
                                        record.samples[0]['AD'], record.samples[0]['DP'])
            snp_dict.update(create_dict)
        elif record.is_indel:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
                                       record.samples[0]['AD'], record.samples[0]['DP'])
            indel_dict.update(create_dict)
        elif record.is_transition:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
                                        record.samples[0]['AD'], record.samples[0]['DP'])
            del_dict.update(create_dict)
        elif record.is_deletion:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
                                        record.samples[0]['AD'], record.samples[0]['DP'])
            transition_dict.update(create_dict)
        else:
            pass
    return snp_dict, indel_dict, del_dict, transition_dict
    

def compute_AF(AD, DP):
    AF = float(AD)/float(DP)
    return float(AF)


def multi_allelic(chrom, pos, ref, alt_list, GT, AD, DP):
    var_dict = defaultdict(list)
    if alt_list[0] !=  None:
        if len(alt_list) == 1:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[genomic_coord].append(value)
        elif len(alt_list) == 2:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[2], DP)
            var_dict[genomic_coord].append(value)
            extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[1])
            extra_value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[extra_coord].append(extra_value)
        elif len(alt_list) == 3:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[3], DP)
            var_dict[genomic_coord].append(value)
            extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[1])
            extra_value = GT, AD,DP, compute_AF(AD[2], DP)
            var_dict[extra_coord].append(extra_value)
            extra_extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[2])
            extra_extra_value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[extra_extra_coord].append(extra_extra_value)
        else:
            pass
    else:
        pass
    return var_dict


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file',
                        help='excel file from workbench export', required=True)
    args = parser.parse_args()
    main(args.vcf_file)
