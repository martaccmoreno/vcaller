from vcaller_funcs.align import *
from vcaller_funcs.processing import *
from vcaller_funcs.variant_calling import *
from vcaller_funcs.compare import *

##################
### MAIN GROUP ###
##################

@click.group()
@click.version_option()
def cli():
    """
    Vcaller, a CLI capable of evoking multiple pre-existing bioinformatics tools, ecletically grouping them into
    commands and subcommands that will perform common variant calling and benchmarking routines.
    """


##########################
### SEQUENCE ALIGNMENT ###
##########################

@cli.group(short_help="Align sequences against reference genome.")
def align():
    """Set of routines to align sequences against a reference genome of choice.
    Read alignment requires two inputs: a reference genome, and reads to align against it.

    In broad strokes, alignment of reads to a reference genome comprises two main steps: Indexing of the reference
    genome (or reads), and the alignment process proper."""


@align.command('bowtie2')
@click.option('--output', '-o', default='bowtie2_output.bam',
              help='Name of the output file (extension will be added automatically)')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bowtie2(output, reference, read1, read2, no_clean):
    """Use the FM-index tool Bowtie 2 for alignment.
    Requires Bowtie 2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""
    if read2 is None:
        func_align_bowtie2(output, reference, read1, no_clean)
    else:
        func_align_bowtie2(output, reference, read1, read2, no_clean)


@align.command('bwa')
@click.option('--output', '-o', default='bwa_out.bam', help='Name of the output file.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, default='')
def align_bwa(output, reference, read1, read2, no_clean):
    """Use the BWA-MEM algorithm for alignment. Requires bwa.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""
    func_align_bwa(output, reference, read1, read2, no_clean)


@align.command('tmap')
@click.option('--output', '-o', default='tmap_output.bam',  # output will have to be redirected
              help='Name of the output file (extension will be added automatically)')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_tmap(output, reference, read1, read2):
    """Use the Ion Torrent-specific aligner TMAP.
    It is only mandatory to include the reference genome file and the reads of a sample as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""
    if read2 is None:
        func_align_tmap(output, reference, read1)
    else:
        func_align_tmap(output, reference, read1, read2)


#######################
### VARIANT CALLING ###
#######################

@cli.group(short_help='Call variants on mapped sequences.')
def call():
    """Set of tools to call variants on files containing sequences previously aligned to a reference genome."""


@call.command('bcftools')
@click.option('--output', '-o', default='bcftools_out.vcf', help='Name of the output file.')
@click.option('--count-orphans', '-A', is_flag=True, help='Count reads with anomolous mate pairs.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_bcftools(output, exome_regions, reference, sample1, sample2, count_orphans, no_clean):
    """Call variants using SAMtools's BCFtools.

    This command calls variants on input aligned sequence files (samples) after calculating their genotype likelihoods.

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """
    func_call_bcftools(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2)


@call.command('freebayes')
@click.option('--output', '-o', default='freebayes_out.vcf',
              help='Name of the output file.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_freebayes(output, exome_regions, reference, sample1, sample2):
    """Call variants using Freebayes.
    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """
    func_call_freebayes(output, exome_regions, reference, sample1, sample2)


@call.command('haplotypecaller')
@click.option('--output', '-o', default='haplotypecaller_out.vcf',
              help='Name of the output file (extension will be added automatically)')
@click.option('--dbsnp', default=None, type=click.Path(exists=True), help='dbSNP file containing a database of '
                                                                          'known SNP IDs.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_haplotypecaller(output, dbsnp, exome_regions, reference, sample1, sample2):
    """Call variants using GATK's HaplotypeCaller.
    The GATK's HaplotypeCaller algorithm is used to call variants on aligned sequence files (samples).

    IMPORTANT NOTES:
    * Sample sequences must have been previously aligned, so that they are in the BAM format.
    * IT IS NECESSARY TO PROCESS THE SAMPLE POST-ALIGNMENT, so that it has read group information.
    *Specifying a dbSNP file with a list of known SNPs is highly recommended.  .

    Through the aforementioned algorithm, this command calls variants on input aligned sequence files (samples),
    simplifying the operation thanks to performing all prerequisite processing steps required by GATK to call variants
    on input files: it ensures that the reference has been indexed through samtools faidx and had a dictionary
    generated through Picard's CreateSequenceDictionary, then it applies the samtools index command on each sample.
    """
    func_call_haplotypecaller(output, dbsnp, exome_regions, reference, sample1, sample2)


@call.command('tvc')
@click.option('--output-dir', '-o', default='.',
              help='Name of output directory; by default save to current directory.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample', type=click.Path(exists=True))
def call_tvc(output_dir, exome_regions, reference, sample):
    """Call variants using TorrentVariantCaller (TVC).
    Only one sample sequence file can be specified. Sample sequence must have been previously aligned, so that it is
    in the SAM/BAM format. A reference genome must be provided.
    """
    func_call_tvc(output_dir, exome_regions, reference, sample)


@call.command('varscan2')
@click.option('--output', '-o', default='bcftools_out.vcf', help='Name of the output file.')
@click.option('--count-orphans', '-A', is_flag=True, help='Count reads with anomolous mate pairs.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_varscan2(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2):
    """Call variants using Varscan2.

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """
    func_call_varscan2(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2)


##### POST-PROCESSING #####
@cli.command('process', short_help='Prepare reads for variant calling.')
@click.option('--output-name', '-o', default=None,
              help='Name of the output file (extension will be added automatically)')
@click.option('--output-dir', '-d', default=None,
              help='Name of output directory; by default save in the same directory as the final output.')
@click.option('--readgroup-info', default='ID:42,PU:123,PL:ILLUMINA,SM:sample1,LB:1',
              type=str, help='Add read group information  to the sample, which MUST follow the format below:\n '
                             'ID:identifier,PU:platform_unit,' '\n'
                             'PL:platform,SM:sample,LB:library' '\n')
@click.option('--add-known-indels', '-i', default='', help='Additional files containing known indel information.',
              multiple=True)
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('known-indels', required=True,
                type=click.Path(exists=True))  # ADD OPTIONS FOR MORE KNOWN INDELS AND SNPS
@click.argument('known-snps', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('sample', required=True, type=click.Path(exists=True))
def process(output_name, output_dir, readgroup_info, add_known_indels, no_clean, known_indels,
            known_snps, reference, sample):
    """Performs a group of steps for the post-processing in preparation for variant calling
    on one SAM/BAM sampl file. A must do for running the gatk subcommand under call."""
    func_process(output_name, output_dir, readgroup_info, add_known_indels, known_indels, known_snps, reference,
                 sample, no_clean)


##### COMPARE #####
@cli.command('compare', short_help='Compare two sets of called variants, '
                                   'with one of them being assumed to be the baseline set.')
@click.option('--output-dir', '-o', default='comparison',  # change to just output
              help='Name of the output directory.')
@click.option('--bed-file', '-b', default=None, help='Only consider variants found in the regions defined by'
                                                     'the provided bed file.')
@click.option('--evaluation-regions', '-e', default=None, help='Define high confidence regions.')
@click.option('--score-field', '-f', default='QUAL',
              help='Choose a custom VCF score field to use as ROC score. Default '
                   'is QUAL.')
@click.option('--sample', '-s', default=None, help='If there is more than 1 sample use to select a sample '
                                                   'or pair of samples.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('baseline', required=True, type=click.Path(exists=True))
@click.argument('calls', required=True, type=click.Path(exists=True))
def compare(output_dir, bed_file, evaluation_regions, score_field, sample, reference, baseline, calls, no_clean):
    """
    Evokes rtgtool's vcfeval to compare a set of baseline calls against a set of query calls, outputting a GA4GH-compliant
    annotated VCF. Next, this VCF is passed into hap.py's qfy.py method in order to compute metrics, namely raw counts
    of TP/FP/FN, as well as their associated precision and recall.
    Because of the way vcfeval's algorithm works, a SDF format file of the reference genome will have to be generated.
    Optionally, regions of the genome wherein to produce the comparison (e.g. exome capture kit regions), as well as
    those considered to be of high confidence, can be defined by providing a BED file.
    To select a pair of samples from each variant data set, use the format <baseline_sample>,<calls_sample>.
    """

    func_compare(output_dir, bed_file, evaluation_regions, score_field, sample,
                 reference, baseline, calls, no_clean)


##### PIPELINE #####
@cli.command('run', short_help='bla')
@click.option('--output', '-o', default='calls.vcf', help='Name of the final VCF output file.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.option('--add-known-indels', '-i', default='', help='Additional files containing known indel information.',
              multiple=True)
@click.option('--readgroup-info', default='ID:42,PU:123,PL:ILLUMINA,SM:sample1,LB:1',
              type=str, help='Add read group information  to the sample, which MUST '
                             'follow the format below:\n'
                             'ID:identifier,PU:platform_unit,' '\n'
                             'PL:platform,SM:sample,LB:library' '\n')
@click.argument('aligner', required=True)
@click.argument('caller', required=True)
@click.argument('known-indels', type=click.Path(exists=True))
@click.argument('known-snps', type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, default='')
def run(output, exome_regions, add_known_indels, readgroup_info, aligner, caller, known_indels, known_snps,
        reference, read1, read2):
    align_out = aligner+'_out.bam'
    process_out = aligner+'_out.processed.bam'

    # align
    if aligner.lower() == 'bowtie2':
        func_align_bowtie2(align_out, reference, read1, read2, no_clean=True)
    elif aligner.lower() == 'bwa':
        func_align_bwa(align_out, reference, read1, read2)
    elif aligner.lower() == 'tmap':
        func_align_tmap(align_out, reference, read1, read2)
    else:
        raise ValueError('The chosen aligner is not valid. Try one of the following: "bowtie2", "bwa", or "tmap".')

    # process
    func_process(process_out, '.', readgroup_info, add_known_indels, known_indels, known_snps, reference, align_out, no_clean=True)

    # call
    if caller.lower() == 'bcftools':
        func_call_bcftools(output, exome_regions, reference, process_out, no_clean=True)
    elif caller.lower() == 'freebayes':
        func_call_freebayes(output, exome_regions, reference, process_out)
    elif caller.lower() == 'haplotypecaller' or caller.lower() == 'gatk':
        func_call_haplotypecaller(output, known_snps, exome_regions, reference, process_out)
    elif caller.lower() == 'tvc':
        func_call_tvc(remove_suffix(output), exome_regions, reference, process_out)
    else:
        raise ValueError('The chosen caller is not valid. Try one of the following: "bcftools", "freebayes", '
                         '"haplotypecaller", or "tvc".')

    return 'Finished!'


