# convert ensembl gff3 to gtf

import numpy as np
import pandas as pd
import csv
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Convert a GFF3 file to basic GTF format')
    parser.add_argument('GFF', metavar='gff3_file', type=str, help='gff3 file to convert')
    parser.add_argument('-o', '--output', type=str, help='file name to write gtf')

    args = parser.parse_args()
    input_file = args.GFF

    # check if output file name specified / sub gff for gtf if no output filename
    if args.output is None:
        output_path = os.path.splitext(input_file)[0] + ".gtf"
    else:
        output_path = args.output

    # check if file already exists
    if os.path.isfile(output_path):
        print(output_path + " already exists")

        overwrite = input("Do you want to overwrite these files? ([Y]/n): ").lower().strip()[:1]
        if not (overwrite == "y" or overwrite == ""):
            sys.exit(1)
        else:
            # remove old files so you don't append new data to old files
            if os.path.isfile(output_path):
                os.remove(output_path)

    # read in gff
    gff3 = pd.read_csv(input_file, sep='\t', comment="#", header=None)
    gff3.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # remove rows that are chromosomes, or CDS/UTR annotations. We don't need these ATM
    gff3 = gff3[~gff3.type.isin(['CDS', 'cds', 'chromosome', 'five_prime_utr', 'three_prime_utr', 'five_prime_UTR', 'three_prime_UTR'])]

    # what does the attributes col start with?
    starts_with = gff3.attributes.str.split(":").str[0]

    # replace any non-standard (i.e. noncoding transcripts) with gene/transcript
    gff3.loc[(starts_with == 'ID=gene'),'type'] = 'gene'
    gff3.loc[gff3.type.str.contains('gene'),'type'] = 'gene'

    gff3.loc[(starts_with == 'ID=transcript'),'type'] = 'transcript'
    gff3.loc[gff3.type.str.contains('RNA'),'type'] = 'transcript'
    gff3.loc[gff3.type.str.contains('exon'),'type'] = 'exon'

    # find gene_id
    gff3['gene_id'] = ''
    gff3.loc[gff3.type=='gene','gene_id'] = gff3.loc[gff3.type=='gene','attributes'].str.split("ID=").str[1].str.split(";").str[0]
    gff3.loc[gff3.type=='transcript','gene_id'] = gff3.loc[gff3.type=='transcript','attributes'].str.split("Parent=").str[1].str.split(";").str[0]

    # find transcript_id
    gff3['transcript_id'] = ''
    gff3.loc[gff3.type=='transcript','transcript_id'] = gff3.loc[gff3.type=='transcript','attributes'].str.split("ID=").str[1].str.split(";").str[0]
    gff3.loc[gff3.type=='exon','transcript_id'] = gff3.loc[gff3.type=='exon','attributes'].str.split("Parent=").str[1].str.split(";").str[0]

    # find exon_id
    gff3['exon_id'] = ''
    # check if ID= or Name= or exon_id=
    replace_with_exon_id = sum(gff3.loc[gff3.type=='exon','attributes'].str.contains('exon_id='))
    replace_with_Name = sum(gff3.loc[gff3.type=='exon','attributes'].str.contains('Name='))
    replace_with_ID = sum(gff3.loc[gff3.type=='exon','attributes'].str.contains('ID='))

    if replace_with_exon_id == len(gff3.loc[gff3.type=='exon','attributes']):
        gff3.loc[gff3.type=='exon','exon_id'] = gff3.loc[gff3.type=='exon','attributes'].str.split("exon_id=").str[1].str.split(";").str[0]
    elif replace_with_Name == len(gff3.loc[gff3.type=='exon','attributes']):
        gff3.loc[gff3.type=='exon','exon_id'] = gff3.loc[gff3.type=='exon','attributes'].str.split("Name=").str[1].str.split(";").str[0]
    elif replace_with_ID == len(gff3.loc[gff3.type=='exon','attributes']):
        gff3.loc[gff3.type=='exon','exon_id'] = gff3.loc[gff3.type=='exon','attributes'].str.split("ID=").str[1].str.split(";").str[0]
    else:
        gff3.loc[gff3.type=='exon','exon_id'] = gff3.loc[gff3.type=='exon','attributes']

    # to find gene_id for exons, copy gene_id from parent transcript
    transcript_ids = gff3.copy()
    transcript_ids = transcript_ids.loc[gff3.type=='transcript',['gene_id','transcript_id']]
    gff3 = pd.merge(gff3, transcript_ids, left_on='transcript_id', right_on='transcript_id', how='left')
    gff3.loc[gff3.type=='exon','gene_id_x'] = gff3.loc[gff3.type=='exon','gene_id_y']

    # add a gtf formatted description/attributes column
    gff3['gtf_desc'] = ''
    gff3.loc[gff3.type=='gene','gtf_desc'] = 'gene_id "' + gff3.loc[gff3.type=='gene','gene_id_x'].astype(str) + '";'
    gff3.loc[gff3.type=='transcript','gtf_desc'] = 'gene_id "' + gff3.loc[gff3.type=='transcript','gene_id_x'].astype(str) + '";' + ' transcript_id "' + gff3.loc[gff3.type=='transcript','transcript_id'].astype(str) + '";'
    gff3.loc[gff3.type=='exon','gtf_desc'] = 'gene_id "' + gff3.loc[gff3.type=='exon','gene_id_x'].astype(str) + '";' + ' transcript_id "' + gff3.loc[gff3.type=='exon','transcript_id'].astype(str) + '";' + ' exon_id "' + gff3.loc[gff3.type=='exon','exon_id'].astype(str) + '";'

    # write file
    gff3.to_csv(output_path, sep='\t', columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase','gtf_desc'], header=False, index=False, quoting=csv.QUOTE_NONE)
