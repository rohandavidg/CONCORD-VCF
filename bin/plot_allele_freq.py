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


#TODO: Dataframe for GATK
def main(dict1, dict2, sample_name1, sample_name2, annotate,logger,
         dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    plot_freq(dict1, dict2,sample_name1, sample_name2, annotate, logger,
              dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_dp_AF(dict1, dict2, sample_name1, sample_name2,annotate, logger,
               dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_dp_bar(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_violin_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                   dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_dist_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                 dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_dist_dp(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                 dict3=False, dict4=False, sample_name3=False, sample_name4=False)
    plot_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
            dict3=False, dict4=False, sample_name3=False, sample_name4=False)
                 

def dict_to_dataframe(some_dict, logger):
    df = pd.DataFrame(some_dict.items(), columns=['position', 'value'])
    df['frequency'] = df.ix[:,1].str[0].str[4]
    df['DP'] = df.ix[:,1].str[0].str[5]
    return df


def dict_to_dataframe_gatk(dict1, dict2, logger):
    fi_df = pd.DataFrame(dict1.items(), columns=['position', 'value'])
    fi_df['frequency'] = fi_df.ix[:,1].str[3]
    fi_df['DP'] = fi_df.ix[:,1].str[2]
    si_df = pd.DataFrame(dict2.items(), columns=['position', 'value'])
    si_df['frequency'] = fi_df.ix[:,1].str[3]
    si_df['DP'] = fi_df.ix[:,1].str[2]
    return fi_df, si_df


def plot_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
            dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, ddict1, logger)
    si_df = dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    fi_df = fi_df.sort_index(by=['frequency'], ascending=[True])
    x1 = np.linspace(0, 1, len(fi_df['frequency']))
    si_df = si_df.sort_index(by=['frequency'],ascending=[True])
    x2 = np.linspace(0, 1, len(si_df['frequency']))
    ax1.scatter(x1, fi_df['frequency'], color='k', marker='x', s=350, vmin=0.,vmax=0.6,label=sample_name1)
    ax1.scatter(x2, si_df['frequency'], color='g', marker='o', s=200, vmin=0.,vmax=0.6,label=sample_name2)
    if dict3:
        ti_df = dict_to_dataframe(dict3, logger)
        x3 = np.linspace(0, 1, len(si_df['ftiquency']))
        ax1.scatter(x1, ti_df['frequency'], color='b', marker='D', s=300, vmin=0.,vmax=0.6,label=sample_name3)        
    if dict4:
        fo_df = dict_to_dataframe(dict4, logger)
        x4 = np.linspace(0, 1, len(fo_df['ftiquency']))
        ax1.scatter(x4, fo_df['frequency'], color='m', marker='d', s=300, vmin=0.,vmax=0.6,label=sample_name4)            
    plt.legend(loc='upper right')
    if annotate:
        annotate_plot(fi_df, ax1, 'total')
        annotate_plot(si_df, ax1, 'total')
    ax1.set_title('Variant frequency correlation', fontsize=20)
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
    plt.savefig('Allele_frequncy_trend.png')

#TODO:dataframe for GATK

def plot_freq(dict1, dict2, sample_name1, sample_name2,annotate, logger,
              dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, logger)
    si_df = dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df.index, fi_df['frequency'], color='r', marker='x', s=350, vmin=0.,vmax=0.6,label=sample_name1)
    ax1.scatter(si_df.index, si_df['frequency'], color='g', marker='o', s=200, vmin=0.,vmax=0.6,label=sample_name2)
    if dict3:
        ti_df = dict_to_dataframe(dict3, logger)
        ax1.scatter(ti_df.index, ti_df['frequency'], color='b', marker='D', s=300, vmin=0.,vmax=0.6,label=sample_name3)
    if dict4:
        fo_df = dict_to_dataframe(dict4, logger)        
        ax1.scatter(fo_df.index, fo_df['frequency'], color='m', marker='d', s=300, vmin=0.,vmax=0.6,label=sample_name4)
    plt.legend(loc='upper left')
    if annotate:
        annotate_plot(fi_df, ax1, 'total')
        annotate_plot(si_df, ax1, 'total')
    ax1.set_title('Variant frequency correlation', fontsize=20)
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('Variant frequency', fontsize=20)
    plt.savefig('Allele_frequency_total_variants.png')


def plot_dp_AF(dict1, dict2, sample_name1, sample_name2,annotate, logger,
               dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, logger)
    si_df =  dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(fi_df['DP'], fi_df['frequency'], color='k', marker='x', s=350, vmin=0., vmax=0.6, label=sample_name1)
    ax1.scatter(si_df['DP'], si_df['frequency'], color='g', marker='o', s=200, vmin=0., vmax=0.6, label=sample_name2)
    if dict3:
        ti_df = dict_to_dataframe(dict3, logger)
        ax1.scatter(ti_df.index, ti_df['frequency'], color='b', marker='D', s=300, vmin=0.,vmax=0.6,label=sample_name3)
    if dict4:
        fo_df = dict_to_dataframe(dict4, logger)
        ax1.scatter(fo_df.index, fo_df['frequency'], color='m', marker='d', s=300, vmin=0.,vmax=0.6,label=sample_name4)
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


def plot_dp_bar(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, logger)
    si_df = dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    ax1 = fig.add_subplot(1,1,1)
    width = 0.50
    reacts1 = ax1.bar(fi_df.index, fi_df['DP'], width, color='green')
    reacts2 = ax1.bar(si_df.index, si_df['DP'], width, color='red')                
    if dict3:
        ti_df = dict_to_dataframe(dict3, logger)
        reacts3 = ax1.bar(ti_df.index, ti_df['DP'], width, color='blue')
    if dict4:
        fo_df = dict_to_dataframe(dict4, logger)
        reacts4 = ax1.bar(fo_df.index, fo_df['DP'], width, color='black')
    plt.legend(loc='upper left')
    ax1.set_title('DP Distribution', fontsize=20)
    ax1.set_xlabel('Total Variants', fontsize=20)
    ax1.set_ylabel('DP', fontsize=20)
    ax1.legend((reacts1[0], reacts2[0]),(sample_name1, sample_name2))
    try:
        plt.savefig('read_depth_distribution.png')
    except ValueError:
        logger.debug("turn off annotation too many variants for plots")

 

def plot_violin_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                   dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1,ddict1, logger)
    si_df = dict_to_dataframe(dict1, ddict2,logger)
    if len(fi_df['frequency']) == len(si_df['frequency']):
        df = pd.DataFrame({sample_name1: fi_df['frequency'], sample_name2: si_df['frequency']})
#        df = pd.concat([fi_df['frequency'], si_df['frequency']], 
        fig = plt.figure(figsize=(20,16))
        sns.set(font_scale=1.8)
        ax = sns.violinplot(df)
        ax.set(xlabel='frequency', ylabel='sample')
        plt.savefig(sample_name1 + "_" + sample_name2 +'_Allele_frequency_distribution.png')
    if dict3:
        ti_df = dict_to_dataframe(dict3, logger)
        if len(fi_df['frequency']) == len(si_df['frequency']) == len(ti_df['frequency']):
            df = pd.DataFrame({sample_name1: fi_df['frequency'], sample_name2: si_df['frequency'],
                               sample_name3: ti_df['frequency']})
            fig = plt.figure(figsize=(20,16))
            sns.set(font_scale=1.8)
            ax = sns.violinplot(df)
            ax.set(xlabel='frequency', ylabel='sample')
            plt.savefig(sample_name1 + "_" + sample_name2 + "_" + sample_name3 + "_" +'_Allele_frequency_distribution.png')
    else:
        pass

    

def plot_dist_af(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                 dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, ddict1,logger)
    si_df = dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    sns.set(font_scale=1.8)
    ax1 = sns.distplot(fi_df.frequency.dropna())
    ax1 = sns.distplot(si_df.frequency.dropna())
    ax1.set(xlabel='frequency', ylabel='sample')
    ax1.set_title(sample_name1 + " vs " + sample_name2)
    plt.legend(loc='upper left')
    if dict3:
        ti_df =  dict_to_dataframe(dict3, logger)
        ax1.set_tile(sample_name1 + " vs " + sample_name2 + " vs " + sample_name3)
        ax1 = sns.distplot(ti_df.frequency.dropna())
    if dict3 and dict4:
        fo_df = dict_to_dataframe(dict3, logger)
        ax1 = sns.distplot(fo_df.frequency.dropna())
    plt.savefig('sampleAllele_frequency_distribution_dist.png')


def plot_dist_dp(dict1, dict2, sample_name1, sample_name2,annotate, logger,
                 dict3=False, dict4=False, sample_name3=False, sample_name4=False):
    fi_df = dict_to_dataframe(dict1, logger)
    si_df = dict_to_dataframe(dict2, logger)
    fig = plt.figure(figsize=(20,16))
    sns.set(font_scale=1.8)
    ax1 = sns.distplot(fi_df['DP'].dropna())
    ax1.set(xlabel='DP', ylabel='sample')
    ax1.set_title("sample " + sample_name1)
    ax2 = sns.distplot(si_df['DP'].dropna())
    ax2.set(xlabel='DP', ylabel='sample')
    ax2.set_title(sample_name1 + " vs " + sample_name2)
    if dict3:
        ti_df =  dict_to_dataframe(dict3, logger)
        ax1.set_tile(sample_name1 + " vs " + sample_name2 + " vs " + sample_name3)
        ax1 = sns.distplot(ti_df['DP'].dropna())
    if dict3 and dict4:
        fo_df = dict_to_dataframe(dict3, logger)
        ax1 = sns.distplot(fo_df['DP'].dropna())
    plt.savefig('sample_DP_dist.png')
    

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
