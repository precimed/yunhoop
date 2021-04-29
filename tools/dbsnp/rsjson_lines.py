import argparse
import json

def printPlacements(info, chrom):
    '''
    rs genomic positions
    '''

    ln=''
    for alleleinfo in info:
        # has top level placement (ptlp) and assembly info
        placement_annot = alleleinfo['placement_annot']
        if alleleinfo['is_ptlp'] and \
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
                    ln="\t".join([chrom, str(pos+1), ref, alt, assembly_name, seq_id, chrom+':'+str(pos+1)])
                    break

    return ln

parser = argparse.ArgumentParser(description='Example of parsing JSON RefSNP Data')
parser.add_argument('-i', dest='input_fn', required=True,
                    help='The name of the input file to parse')
parser.add_argument('-c', dest='chrom', required=True,
                    help='The no. of chromosome')

args = parser.parse_args()

file1 = open(args.input_fn, 'r')
Lines = file1.readlines()
for line in Lines:

    rs_obj = json.loads(line)
    print(rs_obj['refsnp_id'] + "\t")  # rs ID

    if 'primary_snapshot_data' in rs_obj:
        primary_snapshot_data = rs_obj['primary_snapshot_data']
        ln = printPlacements(primary_snapshot_data['placements_with_allele'],args.chrom)
        print('rs'+rs_obj['refsnp_id']+'\t'+ln)
