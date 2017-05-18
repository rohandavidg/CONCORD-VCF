#/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script will asses and  plot the difference is
in allele frequency in the vcf
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


def main(dict1, dict2, logger):
    plot_freq(dict1, dict2, logger)
    plot_dp(dict1, dict2, logger)

def dict_to_dataframe(dict1, dict2, logger):
    fi_df = pd.DataFrame(dict1.items(), columns=['position', 'value'])
    fi_df['frequency'] = fi_df.ix[:,1].str[0].str[4]
    fi_df['DP'] = fi_df.ix[:,1].str[0].str[5]
    si_df = pd.DataFrame(dict2.items(), columns=['position', 'value'])
    si_df['frequency'] = si_df.ix[:,1].str[0].str[4]
    si_df['DP'] = si_df.ix[:,1].str[0].str[5]
    return fi_df, si_df


def plot_freq(dict1, dict2, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df.index, fi_df['frequency'], color='r', marker= 'x', s=100, label='FIRST VCF INPUT')
    ax1.scatter(si_df.index, si_df['frequency'], color='g', marker='o', s=200, label='SECOND VCF INPUT')
    plt.legend(loc='upper left')
    ax1.set_title('Variant frequency correlation', fontsize=20)
    annotate = fi_df['position'].tolist()
    index = fi_df.index.tolist()
    freq = fi_df['frequency'].tolist()
    texts = []
    for i, txt in enumerate(annotate):
        texts.append(ax1.text(index[i],freq[i], txt, rotation=45))
    else:
        pass
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
    plt.savefig('allele_frequency.png')


def plot_dp(dict1, dict2, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df['DP'], fi_df['frequency'], color='r', marker= 'x', s=100, label='FIRST VCF INPUT')
    ax1.scatter(si_df['DP'], si_df['frequency'], color='g', marker='o', s=200, label='SECOND VCF INPUT')
    plt.legend(loc='upper left')
    ax1.set_title('Variant frequency correlation', fontsize=20)
    annotate = fi_df['position'].tolist()
    index = fi_df['DP'].tolist()
    freq = fi_df['frequency'].tolist()
    texts = []
    for i, txt in enumerate(annotate):
        texts.append(ax1.text(index[i],freq[i], txt, rotation=45))
    else:
        pass
    si_annotate = si_df['position'].tolist()
    si_index = si_df['DP'].tolist()
    si_freq = si_df['frequency'].tolist()
    si_texts = []
    for i, txt in enumerate(si_annotate):
        si_texts.append(ax1.text(si_index[i], si_freq[i], txt, rotation=45))
    else:
        pass
    ax1.set_xlabel('DP', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
    plt.savefig('allele_frequency_DP.png')


if __name__ == "__main__":
    main()
