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
import subprocess 

vcf_fixup = "/biotools/biotools/vcflib/2015_3_20/bin/vcffixup"
vcf_remove = "/dlmp/sandbox/cgslIS/rohan/LPEA_CAD_VCF_UTILITIES/VCF_remove.pl"
perl = "/usr/local/biotools/perl/5.10.1/bin/perl"
deannotate = "1000genomes_20130502_GRCh37_nodups.INFO.AC,1000genomes_20130502_GRCh37_nodups.INFO.AF,1000genomes_20130502_GRCh37_nodups.INFO.AFR_AF,1000genomes_20130502_GRCh37_nodups.INFO.AMR_AF,1000genomes_20130502_GRCh37_nodups.INFO.AN,1000genomes_20130502_GRCh37_nodups.INFO.EAS_AF,1000genomes_20130502_GRCh37_nodups.INFO.EUR_AF,1000genomes_20130502_GRCh37_nodups.INFO.SAS_AF,16-CHMKH-A-02-00_GT,17-A3039-A-03-00_GT,17-GL830-A-02-00_GT,AC,ACMG_genes_Dec_2016.name,AF,AN,BlatED.ED,CAVA_ALTANN,CAVA_ALTCLASS,CAVA_ALTSO,CAVA_CLASS,CAVA_CSN,CAVA_ENST,CAVA_GENE,CAVA_IMPACT,CAVA_LOC,CAVA_SO,CAVA_TRINFO,CAVA_TYPE,Clinvar_20160515_GRCh37.all_submitters,Clinvar_20160515_GRCh37.all_traits,Clinvar_20160515_GRCh37.clinical_significance,Clinvar_20160515_GRCh37.conflicted,Clinvar_20160515_GRCh37.pathogenic,Clinvar_20160515_GRCh37.review_status,CompoundHet,CompoundHetEffects,culprit,dbNSFP_v3a_GRCh37.CaddPhred,dbNSFP_v3a_GRCh37.LrtPred,dbNSFP_v3a_GRCh37.MetalrPred,dbNSFP_v3a_GRCh37.MetalrRankscore,dbNSFP_v3a_GRCh37.MetasvmRankscore,dbNSFP_v3a_GRCh37.Polyphen2HdivPred,dbNSFP_v3a_GRCh37.Polyphen2HdivRankscore,dbNSFP_v3a_GRCh37.ProveanConvertedRankscore,dbNSFP_v3a_GRCh37.ProveanPred,dbNSFP_v3a_GRCh37.SiftConvertedRankscore,dbNSFP_v3a_GRCh37.SiftPred,dbNSFP_v3a_GRCh37.Vest3Rankscore,dbSNP139.ID,dbSNP_142_GRCh37p13.ID,DP,ESP_V2_GRCh37.AA._maf,ESP_V2_GRCh37.ALL._maf,ESP_V2_GRCh37.EA._maf,ExAC_r03_GRCh37_nodups.INFO.AC_AFR,ExAC_r03_GRCh37_nodups.INFO.AC_AMR,ExAC_r03_GRCh37_nodups.INFO.AC_EAS,ExAC_r03_GRCh37_nodups.INFO.AC_FIN,ExAC_r03_GRCh37_nodups.INFO.AC_Hemi,ExAC_r03_GRCh37_nodups.INFO.AC_Het,ExAC_r03_GRCh37_nodups.INFO.AC_Hom,ExAC_r03_GRCh37_nodups.INFO.AC_NFE,ExAC_r03_GRCh37_nodups.INFO.AC_SAS,ExAC_r03_GRCh37_nodups.INFO.AF,ExAC_r03_GRCh37_nodups.INFO.AN_NFE,ExAC_r03_GRCh37_nodups.INFO.Hemi_AFR,ExAC_r03_GRCh37_nodups.INFO.Hemi_AMR,ExAC_r03_GRCh37_nodups.INFO.Hemi_EAS,ExAC_r03_GRCh37_nodups.INFO.Hemi_FIN,ExAC_r03_GRCh37_nodups.INFO.Hemi_NFE,ExAC_r03_GRCh37_nodups.INFO.Hemi_OTH,ExAC_r03_GRCh37_nodups.INFO.Hemi_SAS,ExperienceDB.ASSAY,ExperienceDB.CLASS,ExperienceDB.NOTES,ExperienceDB.Submitter,HGMD_2016Q1_GRCh37_nodups.CLASS,HGMD_2016Q1_GRCh37_nodups.PHEN,HGMD_2016Q1_GRCh37_nodups.PubMed,InGeneList,InheritancePattern,LOF,MQ,MQ0,NMD,NS,OMIM_20160427_GRCh37p13_EnsemblGTF_75.Phenotypes,PGx_variants.name,snpeff.Amino_Acid_Change,snpeff.Amino_Acid_length,snpeff.Codon_Change,snpeff.Effect_Impact,snpeff.ERRORS,snpeff.Exon_Rank,snpeff.Functional_Class,snpeff.Gene_Coding,snpeff.Gene_Name,snpeff.Genotype_Number,snpeff.Transcript_BioType,snpeff.Transcript_ID,Link_dbSNP,Link_GoogleSearch,Link_GeneCards,Link_Pubmed,Link_IGV,href"


def main(vcf_file, sample_name, logger, multisample=False, bior=False):
    if bior:
        vcf_fixup_file = convert_vcf_format(vcf_file, logger)
        if multisample:
            vcf_total_dict, vcf_snp_dict, vcf_indel_dict  = vcf_to_dict_multisample(vcf_fixup_file, logger, sample_name)
            return vcf_total_dict, vcf_snp_dict, vcf_indel_dict
    if multisample:
        vcf_total_dict, vcf_snp_dict, vcf_indel_dict  = vcf_to_dict_multisample(vcf_file, logger, sample_name)
        return vcf_total_dict, vcf_snp_dict, vcf_indel_dict
    else:
        vcf_total_dict, vcf_snp_dict, vcf_indel_dict = vcf_to_dict(vcf_file, logger)
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
    

def convert_vcf_format(vcf_file, logger):
    outfile1 = ".".join(vcf_file.split('.')[:-1]) + ".deannotated.vcf"
    outfile2 = ".".join(vcf_file.split('.')[:-1]) + ".deannotated.vcffixup.vcf" 
    cmd = perl + " " + vcf_remove + " -v " + vcf_file + " -o " + deannotate + " > " + outfile1
#    logger.INFO('running cmd {0}'.format(cmd))
    print cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.communicate()
    cmd2 = vcf_fixup + " " + outfile1 + " > " + outfile2
    print cmd2
    p2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out2 = p2.communicate()
    print outfile2
    return outfile2


def vcf_to_dict_multisample(vcf_file, logger, sample_name):
    snp_dict = {}
    indel_dict = {}
    total_dict = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    print sample_name
    print vcf_reader.samples
    sample_index = vcf_reader.samples.index(sample_name)
    for record in vcf_reader:        
        try:
            if record.samples[sample_index]['DP'] and record.samples[sample_index]['AD']:
                key, value = mutli_sample_dict(record.CHROM, record.POS, record.REF, record.ALT, record.samples[sample_index]['GT'],
                                               record.samples[sample_index]['AD'], record.samples[sample_index]['DP'], logger)                
                total_dict[key] = value
                if record.is_snp:
                    key, value = mutli_sample_dict(record.CHROM, record.POS, record.REF, record.ALT, record.samples[sample_index]['GT'],
                                                   record.samples[sample_index]['AD'], record.samples[sample_index]['DP'], logger)
                    snp_dict[key] = value
                elif record.is_indel:
                    key, value = mutli_sample_dict(record.CHROM, record.POS, record.REF, record.ALT, record.samples[sample_index]['GT'],
                                                   record.samples[sample_index]['AD'], record.samples[sample_index]['DP'], logger)
                    indel_dict[key] = value
                else:
                    logger.warning('variant found that is not an indel or snp {0}:{1}{2}>{3}'.format(record.CHROM,  record.POS, record.REF, record.ALT))
        except ValueError:
            pass
    return total_dict, snp_dict, indel_dict


def mutli_sample_dict(chrom, pos, ref, alt_list, GT, AD, DP, logger):
    genomic_coord = str(chrom) + ":g." + str(pos) + str(ref) + ">" + str(alt_list[0])
    value = GT, AD, DP, compute_AF(AD[1], DP)
    return genomic_coord, value


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
#    import argparse
#    from do_logging import configure_logger    
#    parser = argparse.ArgumentParser(description=__doc__)
#    parser.add_argument('-v', dest='vcf_file',
#                        help='excel file from workbench export', required=True)
#    args = parser.parse_args()
#    logger = configure_logger('something')
#    main(args.vcf_file, logger)
    main(vcf_file, sample_name, logger,  multi_sample=False, bior=False)
