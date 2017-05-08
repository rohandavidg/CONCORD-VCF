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
from logging import logging


def main(excel_file, file_name, logging):
    csv_file = parse_excel_template(file_name, excel_file, logging)
    return csv_file


class convert_excel:

    def __init__(self, excel_file, csv_file, index):
        self.index = int(index)
        self.excel_file = excel_file
        self.wb = openpyxl.load_workbook(self.excel_file, read_only=True, data_only=True)
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
                c.writerow([str(cell.value) for cell in r if cell.value and cell.value != '' ])
            except UnicodeEncodeError:
                pass

    def close(self):
        self.csv_file.close()


    @classmethod
    def remove_csv(cls,csv_file):
        return os.remove(csv_file)


def parse_excel_template(file_name, excel_workbook, logging):
    outfile = file_name + ".csv"
    try:
        create_csv_file = convert_excel(excel_workbook, outfile, 0)
        logger.info('creating a csv from excel file with name {0}'.format(file_name))
        new_csv_file = create_csv_file.write_csv()
        create_csv_file.close()
    except:
        logger.debug('failed creating excel file, check excel for errors')
    return outfile


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-e', dest='excel_file',
                        help='excel file from workbench export', required=True)
    parser.add_argument('-n', dest='name', 
                        help='name of the sample or vcf in the excel file', required=True)
    args = parser.parse_args()
    main(args.excel_file, args.name)
