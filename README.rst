========================
how_are_we_stranded_here
========================

.. image:: https://img.shields.io/pypi/v/how_are_we_stranded_here.svg
        :target: https://pypi.python.org/pypi/how_are_we_stranded_here

Python package for testing strandedness of RNA-Seq fastq files


Ever get RNA-Seq data where the library prep or strandedness has been omitted in the methods?

This should save some headaches later in your pipeline and analysis when you realise you've used the wrong strandedness setting (RF/fr-firststrand, FR/fr-secondstrand, unstranded)


Requirements
------------
how_are_we_stranded_here requires the following packages be installed:

kallisto == 0.44.x

python >= 3.6.0

RSeQC

It also requires a transcriptome annotation (.fasta file - e.g. ensembl's .cdna.fasta, or a prebuilt kallisto index), and a corresponding gtf.

Sometimes pseudoalignments will not work with newer versions of kallisto. If this is an issue, we suggest downgrading to 0.44.0.

Installation
------------
.. code-block:: console

        pip install how_are_we_stranded_here


Usage
------------
For basic usage, run check_strandedness with a gtf transcript annotation, transcripts fasta file and fastq read files from one sample.

.. code-block:: console

        check_strandedness --gtf Yeast.gtf --transcripts Yeast_cdna.fasta --reads_1 Sample_A_1.fq.gz --reads_2 Sample_A_2.fq.gz


Output
------------
check_strandedness will print to console the results of infer_experiment.py (http://rseqc.sourceforge.net/#infer-experiment-py), along with an interpretation.

.. code-block:: console

        checking strandedness
        Reading reference gene model stranded_test_WT_yeast_rep1_1_val_1_1/Saccharomyces_cerevisiae.R64-1-1.98.bed ... Done
        Loading SAM/BAM file ...  Total 20000 usable reads were sampled
        This is PairEnd Data
        Fraction of reads failed to determine: 0.0595
        Fraction of reads explained by "1++,1--,2+-,2-+": 0.0073 (0.8% of explainable reads)
        Fraction of reads explained by "1+-,1-+,2++,2--": 0.9332 (99.2% of explainable reads)
        Over 90% of reads explained by "1+-,1-+,2++,2--"
        Data is likely RF/fr-firststrand


Any intermediate files are written to a folder in your current working directory derived from the name of the reads_1 file.


How it Works
------------
check_strandedness.py runs a series of commands to check which direction reads align once mapped in transcripts.

It first creates a kallisto index (or uses a pre-made index) of your organisms transcriptome.

It then maps a small subset of reads (default 200000) to the transcriptome, and uses kallisto's --genomebam argument to project pseudoalignments to genome sorted BAM file.

It finally runs RSeQC's infer_experiment.py to check which direction reads from the first and second pairs are aligned in relation to the transcript strand, and provides output with the likely strandedness of your data.
