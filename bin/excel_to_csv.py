#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This scripts converts an excel to csv
"""

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

log_filename = "CONCORD-VCF"

def main(excel_file_1, sample_name_1, excel_file_2, sample_name_2, annotate_plot):
    logger = configure_logger(log_filename)
    csv_file_1 = parse_excel_template(sample_name_1, excel_file_1, logger)
    csv_file_2 = parse_excel_template(sample_name_2, excel_file_2, logger)
    vcf_dict_1, SNP_dict_1, indel_dict_1 = parse_csv(csv_file_1)
    vcf_dict_2, SNP_dict_2, indel_dict_2 = parse_csv(csv_file_2)
    plot_total = plot_dict.main(vcf_dict_1, vcf_dict_2, "total", logger)
    print annotate_plot
    plot_allele_frequency = plot_allele_freq.main(vcf_dict_1, vcf_dict_2, annotate_plot, logger)
    if SNP_dict_1 and SNP_dict_2:
        plot_SNP = plot_dict.main(SNP_dict_1, SNP_dict_2, "SNP", logger)
    if indel_dict_1 and indel_dict_2:
        plot_indel = plot_dict.main(indel_dict_1, indel_dict_2, "INDEL", logger)

class convert_excel:

    def __init__(self, excel_file, csv_file, index):
        self.index = int(index)
        self.excel_file = excel_file
        print self.excel_file
        self.wb = openpyxl.load_workbook(self.excel_file, read_only=True)
        self.csv_file = open(csv_file,'wb')
        self.writer = csv.writer(self.csv_file, delimiter='\t')


    def __iter__(self):
        return self


    def write_csv(self):
        self.ws = self.wb.worksheets[self.index]
        f = self.csv_file
        c = self.writer
        for r in self.ws.rows:
            try:
                c.writerow([str(cell.value) for cell in r])
            except UnicodeEncodeError:
                pass

    def close(self):
        self.csv_file.close()


    @classmethod
    def remove_csv(cls,csv_file):
        return os.remove(csv_file)


def parse_excel_template(file_name, excel_workbook, logger):
    outfile = file_name + ".csv"
#    try:
    create_csv_file = convert_excel(excel_workbook, outfile, 0)
    logger.info('creating a csv from excel file with name {0}'.format(file_name))
    new_csv_file = create_csv_file.write_csv()
    create_csv_file.close()
#    except:
    logger.debug('failed creating excel file, check excel for errors')
    return outfile


if __name__  ==  "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-FE', dest='excel_file_1',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-FS', dest='sample_name_1', 
                        help='name of the sample or vcf in the excel file', required=True)
    parser.add_argument('-SE', dest='excel_file_2',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-SS', dest='sample_name_2',
                        help='name of the sample or vcf in the excel file', required=True)
    parser.add_argument('-a', dest='annotate_plot', action="store_true",
                        help='flag to turn on or turn off plot annotation with varaint')
    args = parser.parse_args()
    main(args.excel_file_1, args.sample_name_1, args.excel_file_2, 
         args.sample_name_2, args.annotate_plot)
    
