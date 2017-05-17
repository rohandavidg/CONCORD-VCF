#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
this script converts a csv to dict
"""

import pandas as pd
import numpy as np
import re
import pprint
import warnings
warnings.filterwarnings("ignore")

def main(csv_file):
    vcf_dict, snp_dict, indel_dict = parse_csv(csv_file) 
    

def parse_csv(csv_file):
    df = pd.read_csv(csv_file, sep='\t')
    headers = df.columns.values
    keep_headers = headers[6:13]
    new_df = df[keep_headers]
    new_df['chrom'] = new_df.ix[:,1].str.split(':').str[0]
    new_df['genomic'] = new_df.ix[:,1].str.split('.').str[1]
    new_df['pos']  = new_df.ix[:,1].str.split('.').str[1].str.extract('(\d+)').astype(int)
    new_df['ref'] = new_df.ix[:,1].str.split('.').str[1].str.extract('([a-zA-Z]+)').astype(str)
    new_df['alt'] = new_df.ix[:,1].str.split('>').str[1]
    new_df['DP'] = new_df.ix[:,4]
    new_df['AD'] = new_df.ix[:,5]
    new_df['AF'] = new_df.ix[:,6]
    req_headers = list(new_df.columns.values)
    req_headers = ['chrom'] + req_headers[-7:] + ['Exon']
    result_df = new_df[req_headers]
    Total_variants = len(result_df)
    df_indels = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) != 0]
    df_SNP = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) == 0]
    all_dict = create_dict_from_dataframe(result_df)
    snp_dict = create_dict_from_dataframe(df_SNP)
    indel_dict = create_dict_from_dataframe(df_indels)
    return all_dict, snp_dict, indel_dict


def create_dict_from_dataframe(dataframe):
    some_dict = dataframe.groupby('pos')[['chrom', 'ref', 'alt', 'AD', 'AF', 'DP']].apply(lambda x: [x for x in x.values]).to_dict()
    return some_dict


if __name__ == "__main__":
    main(csv_file)
