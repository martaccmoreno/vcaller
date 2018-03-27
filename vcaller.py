import os
from subprocess import run

import click

##### IMPORT CONFIG #####
config = {}
with open(os.path.join(os.path.dirname(__file__), 'config.txt'), 'r') as cfg:
    for line in cfg:
        path = line.split('=')
        config[path[0].strip()] = path[1].strip()


##### AUXILIARY FUNCTIONS #####
def check_existence(filename_list):
    """Check if files with the filenames in the list already exist in the working directory."""
    if sum([os.path.isfile(ifile) for ifile in filename_list]) == len(filename_list):
        return True
    else:
        return False


##### MAIN GROUP #####
@click.group()
@click.version_option()
def cli():
    """Vcaller
    One-stop application that allows the user to leverage several existing bioinformatics tools that constitute variant
    calling pipelines. If desired, the recommended settings for each command can be customized.
    """


##### FILE CHECK COMMANDS #####
@cli.group(short_help="Check file validity.")
def check():
    """Tools to check file integrity and whether or not they can be passed to certain commands."""


@check.command('quickcheck')
@click.argument('file', type=click.Path(exists=True))
def check_quickcheck(file):
    args = ['samtools', 'quickcheck', '-vvv', file]
    run(args)


##### SEQUENCE ALIGNMENT #####

@cli.group(short_help="Align sequences against reference genome.")
def align():
    """Set of tools to align sequences against a reference genome of choice."""


@align.command('bwa')
@click.option('--output', '-o', default='bwa_out.bam', help='Name of the output file.')
@click.option('--nthreads', '-t', default='1', help='Number of CPU threads to use during the alignment step.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(output, nthreads, reference, read1, read2):
    """Use the BWA-MEM algorithm for alignment. Requires bwa.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    # check_index: check if reference genome is already indexed
    suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    if check_existence([reference + suffix for suffix in suffix_list]):
        click.echo('Index files already exist!\n Skipping reference genome indexing.')
    else:
        click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
        run(['bwa', 'index', reference])

    # Align input sequences to the reference genome
    click.echo('Aligning read(s) against the reference genome...')  # find way to make message fancier with custom names
    align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1]
    if read2 is not None:
        align_args += read2
    sam_output = 'bwa_out.sam'
    with open(sam_output, "w") as align_out:
        print(align_args)
        run(align_args, stdout=align_out)

    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', sam_output, '-T', '/tmp/lane_temp', sam_output]
    run(sort_args)

    # clean up intermediary files
    click.echo('Cleaning up %s...' % sam_output)
    run(['rm', sam_output])


@align.command('bowtie2')
@click.option('--output', '-o', default='bowtie2_out.bam',
              help='Name of the output file (extension will be added automatically)')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(output, reference, read1, read2):
    """Use the FM-index tool bowtie2 for alignment. Requires bowtie2.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    ref_basename = reference.split('.')[0]  # bowtie2 uses the reference's basename (no suffix) a lot
    # check_index: check if reference genome is already indexed
    suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    if check_existence([ref_basename + suffix for suffix in suffix_list]):
        click.echo('Index files already exist!\n Skipping reference genome indexing.')
    else:
        click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
        # Index the reference genome
        index_args = [config['bowtie2_path'] + '/bowtie2-build', reference, ref_basename]  # last arg is the
        # output name
        run(index_args)

    # ADD MULTI-THREAD OPTION -p !!
    # Align the input reads against the reference genome
    sam_output = 'bowtie2_out.sam'
    if read2 is None:  # if read is single-ended
        align_args = [config['bowtie2_path'] + '/bowtie2', '-x', ref_basename, read1, '-S', sam_output]
        click.echo('Aligning read %s against the reference genome...' % read1)
    else:  # if paired_end
        align_args = [config['bowtie2_path'] + '/bowtie2', '-x', ref_basename,
                      '-1', read1, '-2', read2, '-S', sam_output]
        click.echo('Aligning reads %s %s against the reference genome...' % (read1, read2))
    run(align_args)

    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', '/tmp/lane_temp', sam_output]
    run(sort_args)

    # clean up intermediary files
    click.echo('Cleaning up %s...' % sam_output)
    run(['rm', sam_output])


##### VARIANT CALLING #####
@cli.group(short_help='Call variants on mapped sequences.')
def call():
    """Set of tools to call variants on files containing sequences previously aligned to a reference genome."""


@call.command('gatk')
@click.option('--output', '-o', default='gatk_out.vcf',
              help='Name of the output file (extension will be added automatically)')
@click.option('--known-snps', '-k', default=None, type=click.Path(exists=True), help='dbSNP file containing a database \
                                                                                     of known SNPs to help improve \
                                                                                     variant calling results.')
@click.option('--clean-up', '-c', default=False, help='Clean up intermediary files to save disk space.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_gatk(output, known_snps, clean_up, reference, sample1, sample2):
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

    # samtools faidx on REFERENCE
    if check_existence([reference + '.fai']):
        click.echo('Reference faidx index file already exists!\nSkipping faidx indexing.')
    else:
        click.echo('Indexing reference file %s...' % reference)
        faidx_args = ['samtools', 'faidx', reference]
        run(faidx_args)

    # generate .dict dictionary file for REFERENCE
    dict_file = reference.split('.')[0] + '.dict'
    if check_existence([dict_file]):
        click.echo(
            'Dictionary file %s already exists.\nSkipping reference genome dictionary file generation.' % dict_file)
    else:
        click.echo('Generating reference genome dictionary %s...' % dict_file)
        dict_vars = ['java', '-jar', config['picard_path'], 'CreateSequenceDictionary', 'R=%s' % reference,
                     'O=%s' % dict_file]
        run(dict_vars)

    sample_list = [sample1] + [s for s in sample2]
    # samtools index on each SAMPLE
    if check_existence([sample + '.bai' for sample in sample_list]):
        click.echo('Sample index .bai files already exist!\nSkipping sample indexing.')
    else:
        click.echo('Need to generate sample index .bai files!\nIndexing sample files %s...' % ', '.join(sample_list))
        index_args = ['samtools', 'index'] + sample_list
        run(index_args)

    # run GATK-HC
    click.echo('Calling variants on samples %s with GATK-HC...' % ', '.join(sample_list))
    if known_snps is None:
        # each sample needs to be preceed by an -I so this is not working for more than one sample?
        gatk_args = [config['gatk4_path'], 'HaplotypeCaller', '-R', reference, '-I'] + sample_list + \
                    ['-O', output]
    else:
        gatk_args = [config['gatk4_path'], 'HaplotypeCaller', '-R', reference, '-I'] + sample_list + \
                    ['--dbsnp', known_snps, '-O', output]
    run(gatk_args)


@call.command('bcftools')
@click.option('--output', '-o', default='bcftools_out.vcf',
              help='Name of the output file (extension will be added automatically)')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_bcftools(output, reference, sample1, sample2):
    """Call variants using SAMtools's BCFtools.

    This command calls variants on input aligned sequence files (samples) after calculating their genotype likelihoods.

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """

    sample_list = [sample1] + [s for s in sample2]
    bcf_output = 'bcftools_out.bcf'
    # bcftools mpileup
    click.echo('Calculating genotype likelihoods for %s...' % ', '.join(sample_list))
    mpileup_args = ['bcftools', 'mpileup', '-Ob', '-o', bcf_output, '-f', reference] + sample_list
    run(mpileup_args)

    # variant calling
    click.echo('Calling variants on %s with BCFtools...' % ', '.join(sample_list))
    call_args = ['bcftools', 'call', '-vmO', 'v', '-o', output, bcf_output]
    run(call_args)

    # clean up intermediary files
    click.echo('Cleaning up %s...' % bcf_output)
    run(['rm', bcf_output])

@call.command('tvc')
@click.option('--output-dir', '-o', default='.',
              help='Name of output directory; by default save to current directory.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_bcftools(output_dir, reference, sample1, sample2):
    """Call variants using TorrentVariantCaller (TVC).

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """

    sample_list = [sample1] + [s for s in sample2]
    click.echo('Calling variants on %s with TVC...' % ', '.join(sample_list))
    call_args = [config['tvc'], '-i', ','.join(sample_list), '-r', reference, '-o', output_dir]
    run(call_args)

##### POST-PROCESSING #####
@cli.command('process', short_help='Prepare reads for variant calling.')
@click.option('--platform', '-p', default='IONTORRENT', help='Platform used in sample sequencing.')
@click.option('--library', '-l', default='library1', help='DNA preparation library identifier.')
@click.option('--sample', '-s', default='NA12878', help='The name of the sample sequenced in a read group.')
# Try to make these 2 options, to see if it's possible to give more than 1 set
@click.argument('known-indels', required=True)
@click.argument('known-snps', required=True)
@click.argument('reference', required=True)
@click.argument('samples', required=True, type=click.Path(exists=True), nargs=-1)
def process(platform, library, sample, known_indels, known_snps, reference, samples):
    """Performs a group of steps for the post-processing in preparation for variant calling
    on one or more SAM/BAM sample files. A must do for running the gatk subcommand under call."""
    sample_list = list(samples)

    # Create a directory structure to store post-processed files, only if it does not exist yet
    # run(['mkdir', '-p', 'post-processing-tmp'])

    ##### STORE MORE STUFF IN THE TMP FOLDER
    # Process each sample at a time
    for smpl in sample_list:
        smpl_name, smpl_extension = smpl.split('.')
        smpl_extension = '.' + smpl_extension

        # sort and convert SAM extension files to BAM
        if smpl_extension.lower() == '.sam':
            click.echo('Sorting and converting %s to BAM...' % smpl)
            sort_args = ['samtools', 'sort', '-O', 'bam', '-o', smpl_name+'.bam', '-T', '/tmp/lane_temp', smpl]
            run(sort_args)
            smpl_extension = '.bam'


        # read groups
        # only works for 1 library atm
        rg_output = smpl_name + '.RG' + smpl_extension
        if not check_existence([rg_output]):
            click.echo('Adding Read Group information to %s...' % smpl)
            read_groups = {'ID': smpl_name, 'PL': platform, 'LB': library, 'PU': 'foo', 'SM': sample}
            rg_args = [config['gatk4_path'], 'AddOrReplaceReadGroups', '-I', smpl,
                       '-O', rg_output, '-RGID', read_groups['ID'], '-RGLB', read_groups['LB'],
                       '-RGPL', read_groups['PL'].upper(), '-RGPU', read_groups['PU'], '-RGSM', read_groups['SM']]
            run(rg_args)
        # NOTE: "ERROR MESSAGE: The platform (ion proton) associated with read group GATKSAMReadGroupRecord
        # @RG:bowtie2_out is not a recognized platform. Allowable options are ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,
        # COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN

        # Mark and remove duplicates (make it an option to not remove, only marking?)
        dup_output = '.'.join(rg_output.split('.')[:-1]) + '.DUP' + smpl_extension
        if not check_existence([dup_output]):
            click.echo('Marking and removing duplicates for %s...' % smpl_name)
            # intermediary file
            dup_args = [config['gatk4_path'], 'MarkDuplicates', '-I', rg_output,
                        '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                        '-M', smpl_name + '.metrics']
            run(dup_args)

        # Realign around indels
        # Using gatk3 because of this
        # https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
        # Preparation: samtools index
        if not check_existence([dup_output + '.bai']):
            click.echo('Indexing %s...' % dup_output)
            run(['samtools', 'index', dup_output])
        # Known indels HAVE to be indexed... how to ensure?
        intervals_output = smpl_name + '.intervals'
        if not check_existence([intervals_output]):
            click.echo('Creating indel realignment intervals for %s...' % smpl_name)
            intervals_args = ['java', '-jar', config['gatk3_path'], '-T', 'RealignerTargetCreator', '-R', reference,
                              '-I', dup_output, '-o', intervals_output, '--known',
                              known_indels]  # only one set of known atm
            run(intervals_args)
        realign_output = '.'.join(dup_output.split('.')[:-1]) + '.RLGN' + smpl_extension
        if not check_existence([realign_output]):
            click.echo('Applying indel realignment based on the intervals for %s...' % smpl_name)
            realign_args = ['java', '-jar', config['gatk3_path'], '-T', 'IndelRealigner', '-R', reference,
                            '-I', dup_output, '-targetIntervals', intervals_output, '-known', known_indels,
                            '-o', realign_output]
            run(realign_args)

        # BQSR
        table_output = smpl_name + '.table'  # should be the ACTUAL name for the file...
        if not check_existence([table_output]):
            click.echo('Creating base score recalibration table for %s...' % smpl_name)
            table_args = [config['gatk4_path'], 'BaseRecalibrator', '-R', reference,
                          '--known-sites', known_snps, '-I', realign_output, '-O', table_output]
            run(table_args)
        bqsr_output = '.'.join(realign_output.split('.')[:-1]) + '.BQSR' + smpl_extension
        if not check_existence([bqsr_output]):
            click.echo('Running base score recalibration on %s...' % smpl_name)
            bqsr_args = [config['gatk4_path'], 'ApplyBQSR',
                         '-I', realign_output, '-bqsr', table_output, '-O', bqsr_output]
            run(bqsr_args)


##### FILTERING #####
@cli.command('filter', short_help='Filter variants in a vcf file.')
@click.argument('reference', required=True)
@click.argument('variants', required=True, type=click.Path(exists=True), nargs=-1)
def process(reference, variants):
    """Desc"""
    variant_list = list(variants)
