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
import seaborn as sns


def main(dict1, dict2, sample_name1, sample_name2, annotate,logger):
    plot_freq(dict1, dict2,sample_name1, sample_name2, annotate, logger)
    plot_dp_AF(dict1, dict2, sample_name1, sample_name2,annotate, logger)
    plot_dp_bar(dict1, dict2, sample_name1, sample_name2,annotate, logger)
    plot_violin_af(dict1, dict2, sample_name1, sample_name2,annotate, logger)
    plot_dist_af(dict1, dict2, sample_name1, sample_name2,annotate, logger)

def dict_to_dataframe(dict1, dict2, logger):
    fi_df = pd.DataFrame(dict1.items(), columns=['position', 'value'])
    fi_df['frequency'] = fi_df.ix[:,1].str[0].str[4]
    fi_df['DP'] = fi_df.ix[:,1].str[0].str[5]
    si_df = pd.DataFrame(dict2.items(), columns=['position', 'value'])
    si_df['frequency'] = si_df.ix[:,1].str[0].str[4]
    si_df['DP'] = si_df.ix[:,1].str[0].str[5]
    return fi_df, si_df


def dict_to_dataframe_gatk(dict1, dict2, logger):
    fi_df = pd.DataFrame(dict1.items(), columns=['position', 'value'])
    fi_df['frequency'] = fi_df.ix[:,1].str[3]
    fi_df['DP'] = fi_df.ix[:,1].str[2]
    si_df = pd.DataFrame(dict2.items(), columns=['position', 'value'])
    si_df['frequency'] = fi_df.ix[:,1].str[3]
    si_df['DP'] = fi_df.ix[:,1].str[2]
    return fi_df, si_df


def plot_freq(dict1, dict2, sample_name1, sample_name2,annotate, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df.index, fi_df['frequency'], color='r', marker='x', s=350, vmin=0.,vmax=0.6,label=sample_name1)
    ax1.scatter(si_df.index, si_df['frequency'], color='g', marker='o', s=200, vmin=0.,vmax=0.6,label=sample_name2)
    plt.legend(loc='upper left')
    if annotate:
        annotate_plot(fi_df, ax1, 'total')
        annotate_plot(si_df, ax1, 'total')
    ax1.set_title('Variant frequency correlation', fontsize=20)
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
#    try:
    plt.savefig('Allele_frequency_total_variants.png')
#    except ValueError:
#        logger.debug("turn off annotation too many variants for plots")
#TODO:dataframe for GATK


def plot_dp_AF(dict1, dict2, sample_name1, sample_name2,annotate, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df['DP'], fi_df['frequency'], color='r', marker='x', s=350, vmin=0., vmax=0.6, label=sample_name1)
    ax1.scatter(si_df['DP'], si_df['frequency'], color='g', marker='o', s=200, vmin=0., vmax=0.6, label=sample_name2)
    plt.legend(loc='upper left')
    if annotate:
        annotate_plot(fi_df, ax1)
        annotate_plot(si_df, ax1)
    ax1.set_title('Variant frequency correlation', fontsize=20)
    ax1.set_xlabel('DP', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
    try:
        plt.savefig('allele_frequency_DP.png')
    except ValueError:
        logger.debug("turn off annotation too many variants for plots")


def plot_dp_bar(dict1, dict2, sample_name1, sample_name2,annotate, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    width = 0.50
    reacts1 = ax1.bar(fi_df.index, fi_df['DP'], width, color='green')
    reacts2 = ax1.bar(si_df.index, si_df['DP'], width, color='red')
    ax1.set_title('DP Distribution', fontsize=20)
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('DP', fontsize=20)
    ax1.legend((reacts1[0], reacts2[0]),(sample_name1, sample_name2))
    try:
        plt.savefig('read_depth_distribution.png')
    except ValueError:
        logger.debug("turn off annotation too many variants for plots")

 

def plot_violin_af(dict1, dict2, sample_name1, sample_name2,annotate, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    if len(fi_df['frequency']) == len(si_df['frequency']):
        fig = plt.figure(figsize=(20,16))
        sns.set(font_scale=1.8)
        ax = sns.violinplot([fi_df['frequency'],si_df['frequency']],horizontal=False)
        ax.set(xlabel='frequency', ylabel='sample')
        plt.savefig(sample_name1 + "_" + sample_name2 +'_Allele_frequency_distribution.png')
    else:
        pass
#    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
#    fig = plt.figure(figsize=(20,16))
#    ax = sns.violinplot(x=fi_df['frequency'])
#    plt.savefig(sample_name1 + '_Allele_frequency_distribution_violin.png')
    


def plot_dist_af(dict1, dict2, sample_name1, sample_name2,annotate, logger):
    fi_df, si_df = dict_to_dataframe(dict1, dict2, logger)
    fig = plt.figure(figsize=(20,16))
    sns.set(font_scale=1.8)
    ax1 = sns.distplot(fi_df.frequency.dropna())
    ax1.set(xlabel='frequency', ylabel='sample')
    ax1.set_title("sample " + sample_name1)
    plt.savefig(sample_name1 + '_Allele_frequency_distribution_dist.png')
    ax2 = sns.distplot(si_df.frequency.dropna())
    ax2.set(xlabel='frequency', ylabel='sample')
    ax2.set_title(sample_name1 + " vs " + sample_name2)
    plt.savefig(sample_name2 + '_Allele_frequency_distribution_dist.png')



def annotate_plot(some_df, plot, total=False):
    annotate = some_df['position'].tolist()
    if total:
        index = some_df.index.tolist()
        freq = some_df['frequency'].tolist()
        texts = []
        for i, txt in enumerate(annotate):
            texts.append(plot.text(index[i], freq[i], txt, rotation=45))
        else:
            pass
    else:
        index = some_df['DP'].tolist()
        freq = some_df['frequency'].tolist()
        texts = []
        for i, txt in enumerate(annotate):
            texts.append(plot.text(index[i], freq[i], txt, rotation=45))
        else:
            pass


if __name__ == "__main__":
    main()
