========================
how_are_we_stranded_here
========================

.. image:: https://img.shields.io/pypi/v/how_are_we_stranded_here.svg
        :target: https://pypi.python.org/pypi/how_are_we_stranded_here


Python package for testing strandedness of RNA-Seq fastq files

Requirements
------------
how_are_we_stranded_here requires the following packages be installed:

kallisto >= 0.44.0
python >= 3.6.0
RSeQC

Installation
------------
.. code-block:: console
        git clone https://github.com/betsig/how_are_we_stranded_here.git
        pip install how_are_we_stranded_here

Usage
------------
For basic usage, run check_strandedness with a gtf transcript annotation, transcripts fasta file and fastq read files from one sample.

.. code-block:: console

        check_strandedness --gtf Yeast.gtf --transcripts Yeast_cdna.fasta --reads_1 Sample_A_1.fq.gz --reads_2 Sample_A_2.fq.gz

How it Works
------------
check_strandedness.py runs a series of commands to check which direction reads align once mapped in transcripts.
It first creates a kallisto index (or uses a pre-made index) of your organisms transcriptome.
It then maps a small subset of reads (default 2 million) to the transcriptome, and uses kallisto's --genomebam argument to project pseudoalignments to genome sorted BAM file.
It finally runs RSeQC's infer_experiment.py to check which direction reads from the first and second pairs are aligned in relation to the transcript strand, and provides output with the likely strandedness of your data.
