
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--gtf", type=str, required=True, help="input gtf file")
    parser.add_argument("--bed", type=str, required=True, help="output bed file")
    parser.add_argument("-t", "--transcript_id_marker", type=str, default="transcript_id", help="text preceeding the transcript id in the 9th field")

    args = parser.parse_args()

    transcript_id_marker = args.transcript_id_marker
    previous_transcript_id = ''

    out_bed = open(args.bed, 'w')

    def write_bed_line():
        starts_array = np.array(starts)
        starts_array = starts_array - 1  # minus one to get the first nt in range
        ends_array = np.array(ends)
        starts_array.sort()
        ends_array.sort()

        block_count = len(starts_array)
        block_sizes = ends_array - starts_array
        block_starts = starts_array - min(starts_array)
        block_sizes = ','.join(str(v) for v in block_sizes.tolist())
        block_starts = ','.join(str(v) for v in block_starts.tolist())

        bed_line = seqname + '\t' + str(min(starts_array)) + '\t' + str(max(ends_array)) + '\t' + previous_transcript_id + '\t' + '0' + '\t' + strand + '\t' + str(min(starts_array)) + '\t' + str(max(ends_array)) + '\t' + '0' + '\t' + str(block_count) + '\t' + block_sizes + ',' + '\t' + block_starts + ','
        out_bed.write(bed_line + '\n')


    # loop through each line, writing bed line  from exon information if a new transcript is encountered.
    #transcript_id_marker = "transcript_id"
    for lines in open(args.gtf):
        # ignore commented lines
        if lines[0] != '#':

            field = lines.strip().split('\t')

            # only use 'exon' lines
            if field[2] == 'exon':
                description = field[8]
                description = description.split(';')

                #clean transcript ids
                transcript_id = [x for x in description if transcript_id_marker in x]
                transcript_id = str(transcript_id[0]).replace(transcript_id_marker, '')
                transcript_id = transcript_id.replace(' ', '')
                transcript_id = transcript_id.replace('\"', '')
                transcript_id = transcript_id.replace('\'', '')

                if transcript_id != previous_transcript_id:

                    # export previous
                    if previous_transcript_id != '':
                        write_bed_line()

                    # start new bed line
                    seqname = field[0]
                    chr_start = str(int(field[3]) - 1)
                    strand = field[6]

                    starts = []
                    ends = []

                starts.append(int(field[3]))
                ends.append(int(field[4]))

                previous_transcript_id = transcript_id

    # write final line
    write_bed_line()
    out_bed.close()
