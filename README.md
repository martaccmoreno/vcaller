Vcaller
========

Variant calling, while simple in principle, can quickly become a daunting task to the
to the average user, as the typical pipeline requires back-and-forth
usage and configuration of multiple command line tools.
Vcaller simplifies this task by serving
as a wrapper which combines multiple pre-existing bioinformatics 
tools into a single package. Vcaller is a command line tool
composed by intuitive building-block commands which retain 
enough flexibility to allow the user to build their own variant 
calling pipeline without having to worry about the minutia 
inherent to the usual "mix-and-match" approach utilized in the field.

What follows is an example pipeline that calls variants on Illumina data:

    $ vcaller align bwa reference.fasta read1.fastq read2.fastq
    $ valler process known_indels.vcf known_snps.vcf reference.fasta sample.bam
    $ vcaller call gatk reference.fasta processed_sample.bam
    
The final output will be a VCF file containing the called variants.
Intermediary step files, while usually cleaned, may optionally be kept.

Vcaller has been tested on a server running Ubuntu 16.04.4 LTS and a laptop
running Linux Mint 18.2 Sonya.

Features
--------

- Be awesome
- Make things faster

Requirements
------------

Software required to run all vcaller features:

* bedtools (https://bedtools.readthedocs.io/en/latest/)
* Bowtie 2*
(http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* GATK3*
(https://software.broadinstitute.org/gatk/download/archive)
* GATK4*
(https://software.broadinstitute.org/gatk/download/)
* Picard*
(https://broadinstitute.github.io/picard/)
* samtools + bcftools
(http://www.htslib.org/download/)
* TVC*
(https://github.com/iontorrent/TS)


Tools marked with a * must have their correct path defined in the config.json file.

Installation
------------

Install vcaller using virtualenv:

    $ cd path/to/vcaller
    $ virtualenv venv
    $ . venv/bin/activate
    $ pip install --editable .
    
Every time you wish to use vcaller command directly from your terminal,
activate the created virtual environment:

    $ . venv/bin/activate
    $ vcaller --help
    
Quickstart
----------

Examples for basic variant calling.

**Aligning to the reference genome:**

*Using bwa*

    $ vcaller align bwa -o bwa_output.bam reference.fa read1.fastq.gz [read2.fastaq.gz]

*Using Bowtie2*

    $ vcaller align bowtie2 -o bt2_output.bam reference.fa read1.fastq.gz [read2.fastq.gz]
    
**Processing aligned reads:**

    $ vcaller process -o processed_output.bam -d path/to/out/dir/ \
    --read-group-info ID:id,PL:SEQUENCINGPLATFORM,PU:platformUnit,SM:sample,LB:library1 \
    -i 1000G_indels.vcf.gz Mills_indels.vcf.gz dbsnp.vcf.gz reference.fa sample.bam

(For more information on Read Group Information, see: 
https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)

(For more information on which known variant databases to use, see:
https://software.broadinstitute.org/gatk/documentation/article.php?id=1247)

**Calling Variants**

*Using GATK4*

    $ vcaller call gatk -o gatk_vars.vcf --dbsnp dbsnp.vcf reference.fa \
    processed_sample.bam

Contribute
----------

- Issue Tracker: github.com/$project/$project/issues
- Source Code: github.com/$project/$project

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: project@google-groups.com

License
-------

The project is licensed under the ??? (BSD?) license.