# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
# Script name: rsjson_demo.py
# Description: a demo script to parse dbSNP RS JSON object.  The script will
# produce tab-delimited output containing tthe assembly version, sequence ID,
# position, reference allele, variant allele and ClinVar clinical significance
# if available.
# Author:  Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------


import os
import argparse
import json
import bz2

parser = argparse.ArgumentParser(description='Example of parsing misc JSON RefSNP Data')
parser.add_argument('-m', dest='merged_fn', required=True,
                    help='Merged snp json file')
parser.add_argument('-n', dest='nosnppos_fn', required=True,
                    help='No snp pos json file')
parser.add_argument('-u', dest='unsupported_fn', required=True,
                    help='Unsupported snp json file')
parser.add_argument('-w', dest='withdrawn_fn', required=True,
                    help='Withdrawn snp json file')
parser.add_argument('-o', dest='other_fn', required=True,
                    help='Other snp json file')
parser.add_argument('-t', dest='out_map_folder', required=True,
                    help='Output map folder')

args = parser.parse_args()

#1) merged refsnps (output: refsnp	merged_into_snp	num_of_merged_into_snps)
if os.path.isfile(args.out_map_folder+'/dbsnp_merge_map.txt'):
    os.remove(args.out_map_folder+'/dbsnp_merge_map.txt')
output_fn = open(args.out_map_folder+'/dbsnp_merge_map.txt', 'a')
input_fn = args.merged_fn
with open(input_fn, 'r') as json_file:
    for line in json_file:
        rs_obj = json.loads(line)

        if 'merged_snapshot_data' in rs_obj:
            merged_into_data = rs_obj['merged_snapshot_data']['merged_into']
            if len(merged_into_data) > 0:
                for i in range(len(merged_into_data)):
                    output_fn.write('rs'+rs_obj['refsnp_id'] + '\t' + 'rs'+merged_into_data[i]+ '\t' + str(len(merged_into_data))+'\n')
output_fn.close()

#2) unmapped/unsupported/withdrawn refsnps (output: refsnp/merged_snp)
tags = ['nosnppos', 'unsupported', 'withdrawn']
for x in tags:
    if os.path.isfile(args.out_map_folder+'/dbsnp_'+x+'.txt'):
        os.remove(args.out_map_folder+'/dbsnp_'+x+'.txt')
    output_fn = open(args.out_map_folder+'/dbsnp_'+x+'.txt', 'a')
    if x == 'nosnppos':
        input_fn = args.nosnppos_fn
    elif (x == 'unsupported'):
        input_fn = args.unsupported_fn
    elif (x == 'withdrawn'):
        input_fn = args.withdrawn_fn

    with open(input_fn, 'r') as json_file:
        for line in json_file:
            rs_obj = json.loads(line)

            output_fn.write('rs'+rs_obj['refsnp_id']+'\n')
            if 'dbsnp1_merges' in rs_obj:
                dbsnp1_merges = rs_obj['dbsnp1_merges']
                if len(dbsnp1_merges) > 0:
                    for i in range(len(dbsnp1_merges)):
                        output_fn.write('rs'+dbsnp1_merges[i]['merged_rsid']+'\n')
    output_fn.close()

#3) other refsnps (output: refsnp	refsnp/merged_snp	refsnp)
if os.path.isfile(args.out_map_folder+'/dbsnp_other_map.txt'):
    os.remove(args.out_map_folder+'/dbsnp_other_map.txt')
output_fn = open(args.out_map_folder+'/dbsnp_other_map.txt', 'a')
input_fn = args.other_fn
with open(input_fn, 'r') as json_file:
    for line in json_file:
        rs_obj = json.loads(line)

        output_fn.write('rs'+rs_obj['refsnp_id']+'\t'+'rs'+rs_obj['refsnp_id']+'\n')
        if 'dbsnp1_merges' in rs_obj:
            dbsnp1_merges = rs_obj['dbsnp1_merges']
            if len(dbsnp1_merges) > 0:
                for i in range(len(dbsnp1_merges)):
                    output_fn.write('rs'+dbsnp1_merges[i]['merged_rsid']+'\t'+'rs'+rs_obj['refsnp_id']+'\n')
    output_fn.close()
