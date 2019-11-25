import numpy as np
import pandas as pd
import csv
import argparse
import os
import sys
import subprocess


parser = argparse.ArgumentParser(description='Check if fastq files are stranded')
parser.add_argument('-g', '--gtf', type=str, help='Genome annotation GTF file', required = True)
parser.add_argument('-fa', '--transcripts', type=str, help='.fasta file with transcript sequences')
parser.add_argument('-n', '--nreads', type=int, help='number of reads to sample', default = 2000000)
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

check_gtf2bed = run_command(cmd = 'python3 gtf2bed.py --help')[1] == b''
if not check_gtf2bed:
    sys.exit("gtf2bed.py is not found in PATH")

check_gff32gtf = run_command(cmd = 'python3 gff32gtf.py --help')[1] == b''
if not check_gff32gtf:
    sys.exit("gff32gtf.py is not found in PATH")

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

#gtf ='/Users/143470/Downloads/Saccharomyces_cerevisiae.R64-1-1.98.gff'

# convert gff to gtf if required
gtf_extension = os.path.splitext(gtf)[1]

if gtf_extension.lower() == '.gff' or gtf_extension.lower() == '.gff3':
    # convert gff to gtf
    gtf_filename = test_folder + '/' + os.path.basename(gtf).replace(gtf_extension, '.gtf')
    cmd = 'python3 gff32gtf.py ' + gtf + ' --output ' + gtf_filename
    print('converting gff to gtf')
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)
else:
    gtf_filename = gtf

# Run gtf2bed
bed_filename = test_folder + '/' + os.path.basename(gtf).replace(gtf_extension, '.bed')
cmd = 'python3 gtf2bed.py --gtf ' + gtf_filename  + ' --bed ' + bed_filename
print('converting gtf to bed')
if print_cmds:
    print('running command: ' + cmd)
subprocess.call(cmd, shell=True)

# make kallisto index
if os.path.exists(kallisto_index_name):
    print('using ' + kallisto_index_name + ' as kallisto index')
else:
    cmd = 'kallisto index -i ' + kallisto_index_name  + ' ' + fasta
    print('generating kallisto index')
    if print_cmds:
        print('running command: ' + cmd)
    subprocess.call(cmd, shell=True)

# take subset of reads
#reads_1 = "/Users/143470/UTS_HPC/SRR1957703_1.fastq.gz"
#reads_2 = "/Users/143470/UTS_HPC/SRR1957703_2.fastq.gz"

print('creating fastq files with first ' + str(n_reads) + ' reads')
reads_1_sample = test_folder + '/' + os.path.basename(reads_1).replace('.fastq' ,'').replace('.fq' ,'').replace('.gz' ,'') + '_sample.fq'
reads_2_sample = test_folder + '/' + os.path.basename(reads_2).replace('.fastq' ,'').replace('.fq' ,'').replace('.gz' ,'') + '_sample.fq'
cmd = 'zcat ' + reads_1 + ' | head -n ' + str(n_reads * 4) + ' > ' + reads_1_sample
if print_cmds:
    print('running command: ' + cmd)
subprocess.call(cmd, shell=True)
if print_cmds:
    print('running command: ' + cmd)
cmd = 'zcat ' + reads_2 + ' | head -n ' + str(n_reads * 4) + ' > ' + reads_2_sample
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

for i in range(4):
    print(result.iloc[i,0])

if float(result.iloc[1,0].replace('Fraction of reads failed to determine: ', '')) > 0.50:
    print('Failed to determine strandedness of > 50% of reads.')
    print('If this is unexpected, try running again with a higher --nreads value')
if float(result.iloc[2,0].replace('Fraction of reads explained by "1++,1--,2+-,2-+": ', '')) > 0.75:
    print('Over 75% of reads explained by "1++,1--,2+-,2-+"')
    print('Data is likely FR/fr-secondstrand')
elif float(result.iloc[3,0].replace('Fraction of reads explained by "1+-,1-+,2++,2--": ', '')) > 0.75:
    print('Over 75% of reads explained by "1+-,1-+,2++,2--"')
    print('Data is likely RF/fr-firststrand')
else:
    print('Less than 75% of reads explained by either category')
    print('Data is likely unstranded')
