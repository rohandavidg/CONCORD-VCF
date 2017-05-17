#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
this script converts a csv to dict
"""

import pandas as pd
import numpy as np
import re
import pprint


def main(csv_file):
    

def parse_csv(csv_file):
    df = pd.read_csv(sep='\t')
    headers = df.columns.values
    keep_headers = headers[6:13]
    new_df = df[keep_headers]
    new_df['chrom'] = new_df['Genomic'].str.split(':').str[0]
    new_df['genomic'] = new_df['Genomic'].str.split('.').str[1]
    new_df['pos']  = new_df['genomic'].str.extract('(\d+)').astype(int)
    new_df['ref'] = new_df['genomic'].str.split('>').str[0].str.extract('([a-zA-Z]+)').astype(str)
    new_df['alt'] = new_df['genomic'].str.split('>').str[1]
    req_aheaders = list(new_df.columns.values)
    req_headers = ['chrom'] + req_headers[-3:] + req_headers[2:7]
    result_df = new_df[req_headers]
    Total_variants = len(result_df)
    df_indels = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) != 0]
    df_SNP = result_df[(result_df.alt.str.len() - result_df.ref.str.len()) == 0]
#TODO:get DP value as well


def create_dict_from_dataframe(dataframe):
    some_dict = dataframe.groupby('pos')[['ref', 'alt', 'Variant Frequency']].apply(lambda x: [x for x in x.values]).to_dict()
    return some_dict
#TODO:get dict

if __name__ == __main__:
    import argparse
    """
    parsing arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='csv_file',
                        help='csv file created from excel export',
                        required=True)

    main(args.csv_file)
