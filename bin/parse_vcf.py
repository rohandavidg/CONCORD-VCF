#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

""" 
parse a vcf file 
tailored for CLC vcf for the moment
"""

import vcf
import os
import csv
from do_logging import configure_logger
from collections import defaultdict
import pprint


def main(vcf_file, logger):
    vcf_total_dict, vcf_snp_dict, vcf_indel_dict  = vcf_to_dict(vcf_file, logger)
    return vcf_total_dict, vcf_snp_dict, vcf_indel_dict


def vcf_to_dict(vcf_file, logger):
    snp_dict = {}
    indel_dict = {}
    total_dict = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
                                   record.samples[0]['AD'], record.samples[0]['DP'], logger)
        total_dict.update(create_dict)
        if record.is_snp:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'], 
                                        record.samples[0]['AD'], record.samples[0]['DP'], logger)
            snp_dict.update(create_dict)
        elif record.is_indel:
            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
                                       record.samples[0]['AD'], record.samples[0]['DP'], logger)
            indel_dict.update(create_dict)
        else:
            logger.warning('variant found that is not an indel or snp {0}:{1}{2}>{3}'.format(record.CHROM,  record.POS, record.REF, record.ALT))
    return total_dict, snp_dict, indel_dict
#        elif record.is_transition:
#            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
#                                        record.samples[0]['AD'], record.samples[0]['DP'], logger)
#            del_dict.update(create_dict)
#        elif record.is_deletion:
#            create_dict = multi_allelic(record.CHROM, record.POS, record.REF, record.ALT, record.samples[0]['GT'],
#                                        record.samples[0]['AD'], record.samples[0]['DP'], logger)
#            transition_dict.update(create_dict)
#        else:
#            pass

    

def compute_AF(AD, DP):
    AF = float(AD)/float(DP)
    return float(AF)


def multi_allelic(chrom, pos, ref, alt_list, GT, AD, DP, logger):
    var_dict = defaultdict(list)
    if alt_list[0] !=  None:
        if len(alt_list) == 1:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[genomic_coord].append(value)
        elif len(alt_list) == 2:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[2], DP)
            var_dict[genomic_coord].append(value)
            extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[1])
            extra_value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[extra_coord].append(extra_value)
        elif len(alt_list) == 3:
            genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
            value = GT, AD, DP, compute_AF(AD[3], DP)
            var_dict[genomic_coord].append(value)
            extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[1])
            extra_value = GT, AD,DP, compute_AF(AD[2], DP)
            var_dict[extra_coord].append(extra_value)
            extra_extra_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[2])
            extra_extra_value = GT, AD, DP, compute_AF(AD[1], DP)
            var_dict[extra_extra_coord].append(extra_extra_value)
        else:
            pass
    else:
        pass
    return var_dict


if __name__ == "__main__":
    main(vcf_file, logger)
