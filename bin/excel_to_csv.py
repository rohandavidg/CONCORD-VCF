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
import pprint
import re
import logging


def main(excel_file):
    


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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-e', dest='excel_file',
                        help='excel file from workbench export', required=True)
    args = parser.parse_args()
    main(args.excel_file)
