#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

""" 
parse a vcf file 
"""

import vcf
import os
import csv
from logging import logging


def main():
    create_vcf_dict = vcf_to_dict(vcf_file)



def vcf_to_dict(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        print record
    





if __name__ == "__main__":
    import argpase
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file',
                        help='excel file from workbench export', required=True)
    args = parser.parse_args()
    main(args.vcf_file)
