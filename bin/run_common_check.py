#/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script is called by find common variants
"""

import vcf
import os
import csv
from do_logging import configure_logger
from collections import defaultdict
import pprint
import parse_vcf
import plot_dict
import plot_allele_freq_vcf
import itertools
import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import scipy.stats as stats
import seaborn as sns
import subprocess
from find_common_variants import *
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s1', dest='sample1',
                        help='create a file with the paths to all the filtered vcf', required=True)
    parser.add_argument('-s2', dest='sample2',
                        help='create a file with the paths to all the filtered vcf', required=True)
    parser.add_argument('-o', dest='outdir',
                        help='outdir to where the output should be located', required=True)
    args = parser.parse_args()
    print "passed main loading analyze_pair"
    analyze_pair(args.sample1,args.sample2,args.outdir)


if __name__ == "__main__":
    main()
