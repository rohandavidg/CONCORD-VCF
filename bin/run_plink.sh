#!/bin/bash


##############################################################
#
#This script runs plink on a set of samples
#
###############################################################

usage ()
{
cat <<EOF
 EOF
##########################################################################################################
##
## Script Options:
##   Required:
##      -D     directory containing vcf
##      -l      enable logging

Usage:

sh run_plink.sh -D $VCF_DIR

EOF
}


PLINK="/usr/local/biotools/plink/1.9/plink"

if [ "$log" == "TRUE"  ]
then
  set -x
fi
