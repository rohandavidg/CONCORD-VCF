#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
This script runs concordance for CLC 
samples
"""

import os
import openpyxl
import csv
import argparse
import collections
from do_logging import configure_logger
from csv_to_dict import parse_csv
import pandas as pd
import plot_dict
import plot_allele_freq


def main(excel_file, sample_name, logger):
    logger = configure_logger(log_filename)
    csv_file = parse_excel_template(sample_name, excel_file, logger)
    return csv_file


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
    create_csv_file = convert_excel(excel_workbook, outfile, 0)
    logger.info('creating a csv from excel file with name {0}'.format(file_name))
    new_csv_file = create_csv_file.write_csv()
    create_csv_file.close()
    logger.debug('failed creating excel file, check excel for errors')
    return outfile


if __name__ == "__main__":
    main(excel_file, sample_name, logger)
