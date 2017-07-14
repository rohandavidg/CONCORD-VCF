#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This scripts converts an excel to csv
"""

import os
import openpyxl
import csv
import argparse
import os
import collections
from do_logging import configure_logger
from csv_to_dict import parse_csv
import pandas as pd
import plot_dict
import plot_allele_freq
import excel_to_csv

log_filename = "CONCORD-VCF"

def main(excel_file_1, sample_name_1,
         excel_file_2, sample_name_2,
         excel_file_3, sample_name_3,
         excel_file_4, sample_name_4,
         annotate_plot):
    logger = configure_logger(log_filename)
    check_samples = def_check_sample_name(sample_name_1, sample_name_2,
                                          sample_name_3,sample_name_4, logger)
    try:
        csv_file_1 = excel_to_csv.main(excel_file_1, sample_name_1, logger)
        csv_file_2 = excel_to_csv.main(excel_file_2,sample_name_2, logger)
        vcf_dict_1, SNP_dict_1, indel_dict_1 = parse_csv(csv_file_1)
        vcf_dict_2, SNP_dict_2, indel_dict_2 = parse_csv(csv_file_2)
        plot_total = plot_dict.main(vcf_dict_1, vcf_dict_2, sample_name_1, sample_name_2, "total", logger)
        plot_allele_frequency = plot_allele_freq.main(vcf_dict_1, vcf_dict_2, sample_name_1, sample_name_2, annotate_plot, logger)
        if SNP_dict_1 and SNP_dict_2:
            plot_snp = plot_dict.main(SNP_dict_1, SNP_dict_2, sample_name_1, sample_name_2, 'SNP', logger)
        if indel_dict_1 and indel_dict_2:
            plot_indel = plot_dict.main(indel_dict_1, indel_dict_2, sample_name_1, sample_name_2, 'INDEL', logger)
        if excel_file_3 and sample_name_3 and excel_file_4 and sample_name_4:
            csv_file_3 = excel_to_csv.main(sample_name_3, excel_file_2, logger)
            csv_file_4 = excel_to_csv.main(sample_name_4, excel_file_2, logger)
            vcf_dict_3, SNP_dict_3, indel_dict_3 = parse_csv(csv_file_3)
            vcf_dict_4, SNP_dict_4, indel_dict_4 = parse_csv(csv_file_4)
            plot_total = plot_dict.main(vcf_dict_1, vcf_dict_2,vcf_dict_3,vcf_dict_4, "total", logger)
            plot_SNP = plot_dict.main(SNP_dict_1, SNP_dict_2, SNP_dict_3, SNP_dict_4, "SNP", logger)
            if indel_dict_3 and indel_dict_4:
                plot_indel = plot_dict.main(indel_dict_1, indel_dict_2, indel_dict_3, indel_dict_4, "INDEL", logger)
    except IndexError:
        print "ERROR: open and save excel sheet as .xlxs"
        sys.exit()

def def_check_sample_name(sample_name_1, sample_name_2, sample_name_3,sample_name_4, logger):
    sample_list = [sample_name_1, sample_name_2, sample_name_3, sample_name_4]
    seen = set([x for x in sample_list if sample_list.count(x) > 1])
    if None in seen:
        pass
    else:
        logger.debug("{0} sample name must be unique".format(seen))



if __name__  ==  "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-FE', dest='excel_file_1',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-FS', dest='sample_name_1', 
                        help='name of the sample or vcf in the excel file', required=True)
    parser.add_argument('-SE', dest='excel_file_2',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-SS', dest='sample_name_2',
                        help='name of the sample or vcf in the excel file', required=True)
    parser.add_argument('-TE', dest='excel_file_3',
                        help='excel file from workbench export')
    parser.add_argument('-TS', dest='sample_name_3',
                        help='name of the sample or vcf in the excel file')
    parser.add_argument('-FOE', dest='excel_file_4',
                        help='excel file from workbench export')
    parser.add_argument('-FOS', dest='sample_name_4',
                        help='name of the sample or vcf in the excel file')    
    parser.add_argument('-a', dest='annotate_plot', action="store_true",
                        help='flag to turn on or turn off plot annotation with varaint')
    args = parser.parse_args()
    main(args.excel_file_1, args.sample_name_1, 
         args.excel_file_2, args.sample_name_2, 
         args.excel_file_3, args.sample_name_3,
         args.excel_file_4, args.sample_name_4,
         args.annotate_plot)
    
