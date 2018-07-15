#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
script converts ngs workbench excel sheet
to a vcf
"""

from collections import OrderedDict
import vcf
import os
import pprint
import argparse
import pandas as pd
import warnings
from csv_to_dict import create_dict_from_dataframe
warnings.filterwarnings("ignore")
import time
import datetime
import subprocess
from operator import itemgetter

today = datetime.date.today()
today_date = str(today).replace("-","")
vcf_fixup = "/biotools/biotools/vcflib/2015_3_20/bin/vcffixup"
vcf_sort = "/usr/local/biotools/vcftools/0.1.8a/bin/vcf-sort"

def main(csv_file, sample_name):
    csv_dict = parse_csv(csv_file)
#    order_dict = OrderedDict(sorted(csv_dict.items(), key=itemgetter([0][0]))
#    order_dict = ordered_dict(csv_dict)
    vcf_file = dict_to_vcf(csv_dict, sample_name)
    fix_vcf = fix_fake_vcf(vcf_file) 

    
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
    all_dict = create_dict_from_dataframe(result_df)
    return all_dict


#def ordered_dict(all_dict):
#    new_dict = {}
#    for k, v in sorted(all_dict.items()):
#        print v[0][0]


def dict_to_vcf(some_csv_dict, sample_name):
    outfile = sample_name + ".vcf"
    with open(outfile, 'w') as fout:
        fout.write("##fileformat=VCFv4.1\n")
        fout.write("##fileDate=" + today_date + "\n")
        fout.write("##fileEncoding=UTF-8\n")
        fout.write("##source=CLC Genomics Grid Worker 6.0 build 60000\n")
        fout.write("##reference=file:/dlmp/clc/data/samples/" + sample_name + '/' + sample_name + "_TCACAGCATT_L001_R1_001%20(Variants).clc\n")
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total number of filtered reads per sample used by variant caller">\n')
        fout.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depth, number of filtered reads supporting the alleles where the first element represents the reference and subsequent elements represent the alternatives in the order listed in the ALT column">\n')
        fout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample_name + " (Variants)\n")
        for k, v in some_csv_dict.items():
            fout.write(str(v[0][0]) + '\t' + str(k) + "\t." + "\t" + str(v[0][1]) +'\t' + v[0][2] + "\t.\t.\t.\t" + "GT:AD:DP\t" + compute_gt(v[0][4]) + ":" + str(get_ref_allele(v[0][3], v[0][5]))+ "," +str(v[0][3]) + ":" +  str(v[0][5]) + '\n')
    return outfile


def fix_fake_vcf(vcf_file):
    new_outfile = vcf_file.split(".")[0] + "_cmb.vcf"
    cmd2 = vcf_fixup + " " + vcf_file + " | " + vcf_sort +  ' -c > ' + new_outfile 
    print cmd2
    p2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out2 = p2.communicate()
    

def get_ref_allele(ref, DP):
    REF_DP = DP - ref
    return REF_DP


def compute_gt(AF):
    if AF < 0.2:
        return "0/0"
    elif 0.2 < AF < 0.8:
        return "0/1"
    elif AF > 0.8:
        return "1/1"
    else:
        print "ERROR: No AF found"



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='csv_file',
                        help='csv_file', required=True)
    parser.add_argument('-s', dest='sample_name',
                        help='sample name of the vcf to create', required=True)
    args = parser.parse_args()
    main(args.csv_file, args.sample_name)
