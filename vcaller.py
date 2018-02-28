import click
from subprocess import run, call, check_output
import os


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

##### PROCESSING #####


##### SEQUENCE ALIGNMENT COMMANDS #####

@cli.group(short_help="Align sequences against reference genome.")
def align():
    """Set of tools to align sequences against a reference genome of choice."""


@align.command('bwa')
@click.option('--name', '-n', default='bwa_out', help='Name of the output file (extension will be added automatically)')
@click.option('--nthreads', '-t', default='1', help='Number of CPU threads to use during the alignment step.')
@click.option('--check_index', '-i', is_flag=True, help='Check if reference is already indexed, and if yes prompt '
                                                        'the user to skip that step.')
@click.option('--sort/--no-sort', default=True, help='Whether the resulting alignment file should be sorted and'
                                                     'converted to BAM.')
# CHECK IF POSSIBLE TO HAVE AN OPTION DEFINED FOR MULTIPLE COMMANDS
@click.option('--clean-up', '-c', default=True, help='Clean up intermediary files to save disk space.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(name, nthreads, check_index, sort, clean_up, reference, read1, read2):
    """Use the BWA-MEM algorithm for alignment. Requires bwa.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    # Index the reference genome; always index by default
    to_index = True

    # (FLAG) check_index: check if already indexed
    if check_index:
        suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        index_files = [reference + suffix for suffix in suffix_list]
        if sum([os.path.isfile(ifile) for ifile in index_files]) == len(suffix_list):
            click.echo('Index files already exist!')
            index_confirmation = input('Are you sure you want to skip the reference genome indexing step? [y/n]: ')
            if index_confirmation.lower() == 'y' or index_confirmation.lower() == 'yes':
                to_index = False
                click.echo('Skipping reference genome indexing.')
        else:
            click.echo('Need to generate index files!')
    if to_index:
        click.echo('Indexing the provided reference genome...')
        run(['bwa', 'index', reference])

    # Align input sequences to the reference genome
    # Missing: Read group info (deal with it in another command?)
    align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1, read2]
    with open(name + '.sam', "w+") as align_out:
        call(align_args, stdout=align_out)
    if sort:
        click.echo('Sorting and converting to BAM...')
        sort_args = ['samtools', 'sort', '-O', 'bam', '-o', name + '.bam', '-T', '/tmp/lane_temp', name + '.sam']
        run(sort_args)

        # clean up
        if clean_up:
            clean_args = ['rm', name + '.sam']
            run(clean_args)

##### VARIANT CALLING #####
@cli.group(short_help='Call variants on mapped sequences.')
def call():
    """Set of tools to call variants on files containing sequences previously aligned to a reference genome."""


@call.command('gatk')
@click.option('--gatk-path', '-g', default=None, type=click.Path(exists=True),
              help='Define the path leading to the GenomeAnalysisTK.jar file. By default, if no path is specified,'
                   'the command will attempt to find the jar file using the "locate" shell command.')
@click.option('--picard-path', '-p', default=None, type=click.Path(exists=True))
@click.option('--name', '-n', default='gatk_out', help='Name of the output file (extension will be added automatically)')
@click.option('--known-snps', '-k', default=None, type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_gatk(gatk_path, picard_path, name, known_snps, reference, sample1, sample2):
    """Use the GATK's HaplotypeCaller algorithm to call variants on aligned sequence files (samples).
    Specifying a dbSNP file with a list of known SNPs is highly recommended.
    Only one sequence file has to be specified. Sequences must have been previously aligned.
    Files must be in the SAM/BAM format. A reference genome must be specified.
    """
    # CHECK IF FILE IS SAM/BAM, ALIGNED, SUCH AND SUCH
    # ...
    # check samtools faidx
    if not os.path.isfile(reference+'.fai'):
        click.echo('Indexing reference file %s...' % reference)
        faidx_args = ['samtools', 'faidx', reference]
        run(faidx_args)

    # check samtools index on each sample (flag to turn it off)
    sample_list = [sample1] + [s for s in sample2]
    sample_indices = [sample+'.bai' for sample in sample_list]
    if not all([os.path.isfile(ifile) for ifile in sample_indices]): # if not all indices exist
        click.echo('Indexing sample input files...')
        index_args = ['samtools', 'index'] + sample_list
    ### Missing .dict file + READGROUPS

    # run GATK-HC
    if gatk_path is None:
        gatk_path = check_output(['locate', 'GenomeAnalysisTK.jar']) # need to check if locate plays nice with "fresh installs"
    gatk_args = ['java', '-jar', gatk_path.strip(), '-R', reference, '-T', 'HaplotypeCaller', '-I'] + sample_list + \
                ['-stand_call_conf', '20', '-o', name+'.vcf']
    click.echo('Ready to ruuuuumble!!!')
    #run(gatk_args)
