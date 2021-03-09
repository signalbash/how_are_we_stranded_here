import numpy as np
import pandas as pd
import csv
import argparse
import os
import sys
import subprocess
import binascii
import re

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'
def main():
    parser = argparse.ArgumentParser(description='Check if fastq files are stranded')
    parser.add_argument('-g', '--gtf', type=str, help='Genome annotation GTF file', required = True)
    parser.add_argument('-fa', '--transcripts', type=str, help='.fasta file with transcript sequences')
    parser.add_argument('-n', '--nreads', type=int, help='number of reads to sample', default = 200000)
    parser.add_argument('-r1', '--reads_1', type=str, help='fastq.gz file (R1)', required = True)
    parser.add_argument('-r2', '--reads_2', type=str, help='fastq.gz file (R2)', required = True)
    parser.add_argument('-k', '--kallisto_index', type=str, help='name of kallisto index (will build under this name if file not found)', default = 'kallisto_index')
    parser.add_argument('-p', '--print_commands', action='store_true', help='Print bash commands as they occur?')

    args = parser.parse_args()
    reads_1 = args.reads_1
    reads_2 = args.reads_2
    n_reads = args.nreads
    kallisto_index_name = args.kallisto_index
    gtf = args.gtf
    fasta = args.transcripts
    print_cmds = args.print_commands

    if fasta is None and (kallisto_index_name is None or not os.path.exists(kallisto_index_name)):
        sys.exit('transcript .fasta sequences are required to generate the kallisto index. Please supply with --transcripts')

    # check if dependancies available
    def run_command(cmd):
        """given shell command, returns communication tuple of stdout and stderr"""
        cmd_result = subprocess.Popen(cmd, shell = True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate()
        return(cmd_result)

    # check that fasta sequence names match bed names
    def check_bed_in_fa(bed_filename, fasta):
        """checks that bed transcript ids match (fuzzy on version numbers) fasta header ids (before space)"""
        cmd = "head " + bed_filename + " | awk '{print $4}' "
        bed_id = run_command(cmd)[0]
        bed_id = bed_id.decode("utf-8").splitlines()

        fa_headers_equal = []
        for i in range(len(bed_id)):
            cmd = 'grep ">" ' + fasta + ' | grep "' + bed_id[i] + '"'
            #print('running command: ' + cmd)
            fa_header = run_command(cmd)[0].decode('utf-8')
            fa_headers_equal.append(re.sub("[.][0-9]*", "", fa_header.split(" ")[0].replace(">", "")) == bed_id[i] or fa_header.split(" ")[0].replace(">", "") == bed_id[i])

        return sum(fa_headers_equal) > 0


    check_gtf2bed = run_command(cmd = 'gtf2bed --help')[1] == b''
    if not check_gtf2bed:
        sys.exit("gtf2bed is not found in PATH")

    check_gff32gtf = run_command(cmd = 'gff32gtf --help')[1] == b''
    if not check_gff32gtf:
        sys.exit("gff32gtf is not found in PATH")

    check_kallisto = run_command(cmd = 'kallisto version')
    if not check_kallisto[1] == b'':
        sys.exit("kallisto is not found in PATH. Please install from https://pachterlab.github.io/kallisto")
    else:
        kallisto_version = str(check_kallisto[0]).split('version ')[1].replace("\\n'","")
        if int(kallisto_version.split('.')[1]) < 44:
            sys.exit('Found kallisto ' + kallisto_version +' , but version >= 0.44.0 is required. Please install from https://pachterlab.github.io/kallisto')

    check_RSeQC = run_command(cmd = 'infer_experiment.py --help')[1] == b''
    if not check_RSeQC:
        sys.exit("infer_experiment.py (RSeQC) is not found in PATH. Please install from http://rseqc.sourceforge.net/#installation")


    # make a test_folder
    test_folder = 'stranded_test_' + os.path.basename(reads_1).replace('.fastq' ,'').replace('.fq' ,'').replace('.gz' ,'')
    if not os.path.isdir(test_folder):
        # make directory
        os.mkdir(test_folder)
    else:
        index_n = 1
        while os.path.isdir(test_folder + '_' + str(index_n)):
            index_n = index_n + 1
        #make directory
        test_folder = test_folder + '_' + str(index_n)
        os.mkdir(test_folder)
    print("Results stored in: " + test_folder)

    # convert gff to gtf if required
    gtf_extension = os.path.splitext(gtf)[1]

    if gtf_extension.lower() == '.gff' or gtf_extension.lower() == '.gff3':
        # convert gff to gtf
        gtf_filename = test_folder + '/' + os.path.basename(gtf).replace(gtf_extension, '.gtf')
        cmd = 'gff32gtf ' + gtf + ' --output ' + gtf_filename
        print('converting gff to gtf')
        if print_cmds:
            print('running command: ' + cmd)
        subprocess.call(cmd, shell=True)
    else:
        gtf_filename = gtf

    # Run gtf2bed

    bed_filename = test_folder + '/' + os.path.basename(gtf).replace(gtf_extension, '.bed')
    cmd = 'gtf2bed --gtf ' + gtf_filename  + ' --bed ' + bed_filename
    print('converting gtf to bed')
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)

    # make kallisto index
    if os.path.exists(kallisto_index_name):
        print('using ' + kallisto_index_name + ' as kallisto index')
    else:
        print('Checking if fasta headers and bed file transcript_ids match...')
        check_bed = check_bed_in_fa(bed_filename, fasta)
        if not check_bed:
            print("Can't find transcript ids from " + fasta + " in " + bed_filename)
            print("Trying to converting fasta header format to match transcript ids to the BED file...")
            cmd = "sed 's/[|]/ /g' " + fasta + " > " + test_folder + "/transcripts.fa"
            if print_cmds:
                print('running command: ' + cmd)
            subprocess.call(cmd, shell=True)
            fasta = test_folder + "/transcripts.fa"
            check_bed_converted = check_bed_in_fa(bed_filename, fasta)
            if not check_bed_converted:
                subprocess.call("rm -f " + fasta, shell=True)
                sys.exit("Can't find any of the first 10 BED transcript_ids in fasta file... Check that these match")
        else:
            print("OK!")

        cmd = 'kallisto index -i ' + kallisto_index_name  + ' ' + fasta
        print('generating kallisto index')
        if print_cmds:
            print('running command: ' + cmd)
        subprocess.call(cmd, shell=True)

        if not check_bed:
            cmd = "rm -f " + fasta
            subprocess.call(cmd, shell=True)

    print('creating fastq files with first ' + str(n_reads) + ' reads')
    reads_1_sample = test_folder + '/' + os.path.basename(reads_1).replace('.fastq' ,'').replace('.fq' ,'').replace('.gz' ,'') + '_sample.fq'
    reads_2_sample = test_folder + '/' + os.path.basename(reads_2).replace('.fastq' ,'').replace('.fq' ,'').replace('.gz' ,'') + '_sample.fq'
    # check if the fasta is gzipped
    if(is_gz_file(reads_1)):
        cmd = 'zcat < ' + reads_1 + ' | head -n ' + str(n_reads * 4) + ' > ' + reads_1_sample
    else:
        cmd = 'head ' + reads_1 + ' -n ' + str(n_reads * 4) + ' > ' + reads_1_sample
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)
    if print_cmds:
        print('running command: ' + cmd)
    # check if the fasta is gzipped
    if(is_gz_file(reads_2)):
        cmd = 'zcat < ' + reads_2 + ' | head -n ' + str(n_reads * 4) + ' > ' + reads_2_sample
    else:
        cmd = 'head ' + reads_2 + ' -n ' + str(n_reads * 4) + ' > ' + reads_2_sample
    subprocess.call(cmd, shell=True)

    # align with kallisto
    print('quantifying with kallisto')
    cmd = 'kallisto quant -i  ' + kallisto_index_name + '  -o ' + test_folder + '/' + 'kallisto_strand_test --genomebam --gtf ' + gtf_filename + ' ' + reads_1_sample + ' ' + reads_2_sample
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)

    # check strandedness w/ 2million alignments
    #n_reads = 2000000
    print('checking strandedness')
    cmd = 'infer_experiment.py -r ' + bed_filename + ' -s ' + str(n_reads) + ' -i ' + test_folder + '/' + 'kallisto_strand_test/pseudoalignments.bam > ' + test_folder + '/' + 'strandedness_check.txt'
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)

    result = pd.read_csv(test_folder + '/' + 'strandedness_check.txt', sep="\n", header=None)

    failed = float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', ''))
    fwd = float(result.iloc[2,0].replace('Fraction of reads explained by "1++,1--,2+-,2-+": ', ''))
    rev = float(result.iloc[3,0].replace('Fraction of reads explained by "1+-,1-+,2++,2--": ', ''))
    fwd_percent = fwd/(fwd+rev)
    rev_percent = rev/(fwd+rev)

    print(result.iloc[0,0])
    print(result.iloc[1,0])
    print(result.iloc[2,0] + " (" + str(round(fwd_percent*100, 1)) + "% of explainable reads)")
    print(result.iloc[3,0] + " (" + str(round(rev_percent*100, 1)) + "% of explainable reads)")


    if float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', '')) > 0.50:
        print('Failed to determine strandedness of > 50% of reads.')
        print('If this is unexpected, try running again with a higher --nreads value')
    if fwd_percent > 0.9:
        print('Over 90% of reads explained by "1++,1--,2+-,2-+"')
        print('Data is likely FR/fr-secondstrand')
    elif rev_percent > 0.9:
        print('Over 90% of reads explained by "1+-,1-+,2++,2--"')
        print('Data is likely RF/fr-firststrand')
    elif max(fwd_percent, rev_percent) < 0.6:
        print('Under 60% of reads explained by one direction')
        print('Data is likely unstranded')
    else:
        print('Data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6)')
        print('Please check your data for low quality and contaminating reads before proceeding')
