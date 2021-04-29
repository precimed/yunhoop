import argparse
import json

def printPlacements(info, chrom, assemb):
    '''
    rs genomic positions
    '''

    ln=''
    for alleleinfo in info:
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
                        break
    return ln

parser = argparse.ArgumentParser(description='Example of parsing JSON RefSNP Data')
parser.add_argument('-i', dest='input_fn', required=True,
                    help='The name of the input file to parse')
parser.add_argument('-c', dest='chrom', required=True,
                    help='The no. of chromosome')

args = parser.parse_args()

with open(args.input_fn) as f:
    rs_obj = json.load(f)

    if 'primary_snapshot_data' in rs_obj:
        primary_snapshot_data = rs_obj['primary_snapshot_data']
        ln = printPlacements(primary_snapshot_data['placements_with_allele'],args.chrom, 'GRCh38')
        print('rs'+rs_obj['refsnp_id']+'\t'+ln+'\t'+primary_snapshot_data['variant_type'].upper())
        ln = printPlacements(primary_snapshot_data['placements_with_allele'],args.chrom, 'GRCh37')
        print('rs'+rs_obj['refsnp_id']+'\t'+ln+'\t'+primary_snapshot_data['variant_type'].upper())
