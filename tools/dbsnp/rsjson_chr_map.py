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


def printAllele_annotations(primary_refsnp):
    '''
    rs clinical significance
    '''
    for annot in primary_refsnp['allele_annotations']:
        for clininfo in annot['clinical']:
            print(",".join(clininfo['clinical_significances']))

def printPlacements(info, chrom, assemb):
    '''
    rs genomic positions
    '''

    ln=''
    for alleleinfo in info:
        if 'NC_' not in alleleinfo['seq_id']: continue
        # has top level placement (ptlp) and assembly info
        placement_annot = alleleinfo['placement_annot']
        #if alleleinfo['is_ptlp'] and \
        if placement_annot['seq_type'] == 'refseq_chromosome' and \
                len(placement_annot['seq_id_traits_by_assembly']) > 0:
            assembly_name = placement_annot[
                'seq_id_traits_by_assembly'][0]['assembly_name']

            for a in alleleinfo['alleles']:
                spdi = a['allele']['spdi']
                if spdi['inserted_sequence'] != spdi['deleted_sequence']:
                    (ref, alt, pos, seq_id) = (spdi['deleted_sequence'],
                                               spdi['inserted_sequence'],
                                               spdi['position'],
                                               spdi['seq_id'])
                    if assemb in assembly_name:
                        ln = "\t".join([chrom+':'+str(pos+1), chrom, str(pos), ref, alt, assembly_name, seq_id, chrom+':'+str(pos)])
                        return ln
    return ln

parser = argparse.ArgumentParser(description='Example of parsing JSON RefSNP Data')
parser.add_argument('-i', dest='input_fn', required=True,
                    help='The name of the input file to parse')
parser.add_argument('-c', dest='chrom', required=True,
                    help='The no. of chromosome')
parser.add_argument('-o', dest='output_pf', required=True,
                    help='The prefix of the output file')

args = parser.parse_args()

if os.path.isfile(args.output_pf+"_b38.txt"):
    os.remove(args.output_pf+"_b38.txt")
if os.path.isfile(args.output_pf+"_b37.txt"):
    os.remove(args.output_pf+"_b37.txt")
output_fn1 = open(args.output_pf+"_b38.txt", 'a')
output_fn2 = open(args.output_pf+"_b37.txt", 'a')
with bz2.BZ2File(args.input_fn, 'rb') as f_in:
    for line in f_in:
        rs_obj = json.loads(line.decode('utf-8'))
        #print(rs_obj['refsnp_id'] + "\t")  # rs ID
        #json_formatted_str = json.dumps(rs_obj, indent=2)
        #print(json_formatted_str)

        if 'primary_snapshot_data' in rs_obj:
            primary_snapshot_data = rs_obj['primary_snapshot_data']
            ln = printPlacements(primary_snapshot_data['placements_with_allele'],args.chrom, 'GRCh38')
            if ln != '':
                output_fn1.write('rs'+rs_obj['refsnp_id']+'\t'+ln+'\t'+primary_snapshot_data['variant_type'].upper()+'\t'+primary_snapshot_data['anchor']+'\n')
            ln = printPlacements(primary_snapshot_data['placements_with_allele'],args.chrom, 'GRCh37')
            if ln != '':
                output_fn2.write('rs'+rs_obj['refsnp_id']+'\t'+ln+'\t'+primary_snapshot_data['variant_type'].upper()+'\t'+primary_snapshot_data['anchor']+'\n')

output_fn1.close()
output_fn2.close()
