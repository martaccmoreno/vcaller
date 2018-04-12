Vcaller
========

Variant calling, while simple in principle, can quickly become a daunting task to the
to the average user, as the typical pipeline requires back-and-forth
usage and configuration of multiple command line tools.
Vcaller (provisory name) seeks to simplify this task by serving
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

Features
--------

- Be awesome
- Make things faster

Installation
------------

Install vcaller using virtualenv:

    $ cd path/to/vcaller
    $ virtualenv venv
    $ . venv/bin/activate
    $ pip install --editable .

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