#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script compares the concordance
between a vcf and workbench excel
"""

import vcf
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
import numpy as np
import re
import pprint
import warnings
from  parse_vcf import multi_allelic
warnings.filterwarnings("ignore")
import random 


def main(excel_file, vcf_files, test_name):
    random_number = random.randint(0, 10000) 
    log_filename = "CONCORD-VCF-EXCEL_" + test_name + "_" + str(random_number)
    logger = configure_logger(log_filename)
    csv_file = excel_to_csv.main(excel_file, test_name+'_' +str(random_number), logger)
    variant_dict = parse_csv(csv_file)
    test_vcf_dict = parse_list_vcf(vcf_files, logger)
    generate_results = compare_dicts(variant_dict, test_vcf_dict, test_name, random_number,logger)


def parse_csv(csv_file):
    df = pd.read_csv(csv_file, sep='\t')
    headers = df.columns.values
    keep_headers = headers[1:13]
    new_df = df[keep_headers]
    new_df['Sample'] = new_df.ix[:,0]
    new_df['chrom'] = new_df.ix[:,6].str.split(':').str[0]
    new_df['genomic'] = new_df.ix[:,6].str.split('.').str[1]
    new_df['pos']  = new_df.ix[:,6].str.split('.').str[1].str.extract('(\d+)').astype(int)
    new_df['ref'] = new_df.ix[:,6].str.split('.').str[1].str.extract('([a-zA-Z]+)').astype(str)
    new_df['alt'] = new_df.ix[:,6].str.split('>').str[1]
    new_df['DP'] = new_df.ix[:,9]
    new_df['AD'] = new_df.ix[:,10]
    new_df['AF'] = new_df.ix[:,11]
    new_df['key'] = new_df['pos'].apply(str) + "_" +  new_df['Sample'] + "_" + new_df['ref'] + "_" +  new_df['alt']
    req_headers = list(new_df.columns.values)
    req_headers = ['chrom'] + req_headers[-7:] + ['Exon']
    result_df = new_df[req_headers]
    Total_variants = len(result_df)
    df_indels = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) != 0]
    df_SNP = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) == 0]
    all_dict = create_dict_from_dataframe(result_df)
    snp_dict = create_dict_from_dataframe(df_SNP)
    indel_dict = create_dict_from_dataframe(df_indels)
    return all_dict


def create_dict_from_dataframe(dataframe):
    some_dict = dataframe.groupby('key')[['chrom', 'pos', 'ref', 'alt', 'AD', 'AF', 'DP']].apply(lambda x: [x for x in x.values]).to_dict()
    return some_dict


def parse_list_vcf(vcf_files, logger):
    vcf_dict = {}
    for vcf in vcf_files:
        logger.info("looking into vcf {0}".format(vcf))
        vcf_sample = get_sample_name(vcf)
        vcf_run_name = get_run_name(vcf) 
        sample_dict = vcf_to_dict_siglesample(vcf, logger, vcf_sample, vcf_run_name)
        vcf_dict.update(sample_dict)
    return vcf_dict

def vcf_to_dict_siglesample(vcf_file, logger, sample_name, vcf_run_name):
    total_dict = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        key = str(record.POS) + "_" + sample_name + "_" + record.REF + "_" + str(record.ALT[0])
        try:
            value = [record.CHROM + ':' + str(record.POS) +  record.REF + '>' + str(record.ALT[0]), vcf_run_name, record.INFO['CAVA_CSN']]
            total_dict[key] = value
        except:
            value = [record.CHROM + ':' + str(record.POS) +  record.REF + '>' + str(record.ALT[0]), vcf_run_name]
            total_dict[key] = value
    return total_dict


def get_sample_name(vcf):
    sample_name =  vcf.split('/')[-2]
    return sample_name

def get_run_name(vcf):
    run_name =  vcf.split('/')[-4]
    return run_name


def compare_dicts(all_dict, vcf_dict, test, random, logger):
    with open("missing_NGSworkbench_variant_" + test + str(random) + ".txt", 'w') as fout:
        for key, value in vcf_dict.items():
            try:
                all_dict[key]
            except KeyError:
                try:
                    out = [value[0], test, key.split("_")[1], value[1], value[2]]
                    fout.write("\t".join(str(i) for i in out) + "\n")
                    logger.debug("{0} is missing in {1} for sample {2}, found in vcf but not in the excel".format(value,test, 
                                                                                                                  key.split("_")[1]))
                except IndexError:
                    out = [value[0], test, key.split("_")[1], value[1]]
                    fout.write("\t".join(str(i) for i in out) + "\n")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-e', dest='excel_file',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-v', dest='vcf_files', nargs="+", 
                        help='list of vcf file from that run', required=True)
    parser.add_argument('-t', dest='test_name',
                        help='name of the panel', required= True)
    args = parser.parse_args()
    main(args.excel_file, args.vcf_files, args.test_name)
