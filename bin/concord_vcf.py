#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
check concordance between VCF
"""

import vcf
import os
import csv
from do_logging import configure_logger
from collections import defaultdict
import pprint
import parse_vcf
import plot_dict
import plot_allele_freq


log_filename = "CONCORD-VCF"


def main(VCF_file_1, sample_name_1,
         VCF_file_2, sample_name_2, annotate_plot):
    logger = configure_logger(log_filename)
    VCF_file1_total, VCF_file1_SNP_dict, VCF_file1_indel_dict = parse_vcf.main(VCF_file_1, logger)
    VCF_file2_total, VCF_file2_SNP_dict, VCF_file2_indel_dict = parse_vcf.main(VCF_file_2, logger)
    plot_total = plot_dict.main(VCF_file1_total, VCF_file2_total, "total", logger)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-FV', dest='VCF_file_1',
                        help='VCF file from workbench export', required=True)
    parser.add_argument('-FS', dest='sample_name_1',
                        help='name of the sample or vcf in the VCF file', required=True)
    parser.add_argument('-SV', dest='VCF_file_2',
                        help='VCF file from workbench export', required=True)
    parser.add_argument('-SF', dest='sample_name_2',
                        help='name of the sample or vcf in the VCF file', required=True)
    parser.add_argument('-a', dest='annotate_plot', action="store_true",
                        help='flag to turn on or turn off plot annotation with varaint')
    args = parser.parse_args()
    main(args.VCF_file_1, args.sample_name_1,
         args.VCF_file_2, args.sample_name_2, args.annotate_plot)

