#/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
filters vcf based on AF, AD and ROI
"""

import argparse
import collections
import os
import vcf
import pprint
import time
import datetime
import subprocess
import logging

INTERSECTBED = "/usr/local/biotools/bedtools/2.20.1/bin/intersectBed"
VCF_FILTER = "/dlmp/sandbox/cgslIS/rohan/LPEA_CAD_VCF_UTILITIES/VCF_filter.pl"
PERL = "/usr/local/biotools/perl/5.10.1/bin/perl"
filename = "filter_vcf"


def main():
    args = parse_args()
    run(args.vcf_file, args.bed_file, args.allele_frequency, 
        args.allele_depth, args.sample_name)


def run(vcf_file, bed_file, allele_frequency,
        allele_depth, sample_name):
    header = create_tmp_header(vcf_file, sample_name)
    intersect_tmpfile = intesect_vcf_bed(vcf_file, bed_file, sample_name)
    ROI_filterd_vcf = cat_files(intersect_tmpfile, sample_name)
    filter_vcf_thresholds = vcf_filter_filtering(ROI_filterd_vcf, allele_depth, allele_frequency,
                                                 sample_name)
    
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file',
                        help='clc vcf file', required=True)
    parser.add_argument('-b', dest='bed_file',
                        help='bed file of the test', required=True)
    parser.add_argument('-af', dest='allele_frequency',
                        help='allele frequency for the variant to be filtered on', 
                        required=True)
    parser.add_argument('-ad', dest='allele_depth',
                        help='allele depth for the variant to be filtered on', 
                        required=True)
    parser.add_argument('-s', dest='sample_name',
                        help='sample name of vcf', required=True)
    args = parser.parse_args()
    return args


def configure_logger(filename):
    """
    setting up logging
    """
    logger = logging.getLogger(filename)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime(filename+"-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def create_tmp_header(vcf_file, sample_name):
    with open(vcf_file) as fin, open(sample_name + "_header.vcf", 'w') as fout:
        vcf_line = fin.readlines()
        for raw_line in vcf_line:
            line = raw_line.strip()
            if line.startswith("#"):
                fout.write(line + '\n')


def intesect_vcf_bed(vcf_file, bed_file, sample_name):
    tmp_out =  sample_name + "_ROI_filtered.tmp"
    cmd = INTERSECTBED + ' -a ' + vcf_file + ' -b ' + bed_file + " > " + tmp_out
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out = p.communicate()
    return tmp_out


def cat_files(tmp_out, sample_name):
    vcf_header = sample_name + "_header.vcf"
    outfile = sample_name + "_ROI_filtered.vcf"
    with open(outfile, 'w') as fout:
        with open(vcf_header) as hin, open(tmp_out) as tin:
            fout.write(hin.read())
            fout.write(tin.read())
    return outfile


def vcf_filter_filtering(outfile, ad, af, sample_name):
    filtered_vcf = sample_name + "_filtered.vcf"
    cmd = PERL + " " + VCF_FILTER + " -v " + outfile + " -a " + ad + " -f " + af + " > " + filtered_vcf
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out = p.communicate()


if __name__ == "__main__":
    main()
