import click
from subprocess import run, call
import os


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
@click.option('--name', '-n', default='bwa_out', help='Name of the output file (extension will be added automatically)')
@click.option('--nthreads', '-t', default='1', help='Number of CPU threads to use during the alignment step.')
@click.option('--check_index', '-i', is_flag=True, help='Check if reference is already indexed, and if yes prompt '
                                                        'the user to skip that step.')
# CHECK IF POSSIBLE TO HAVE AN OPTION DEFINED FOR MULTIPLE COMMANDS
@click.option('--clean-up', '-c', default=False, help='Clean up intermediary files to save disk space.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(name, nthreads, check_index, clean_up, reference, read1, read2):
    """Use the BWA-MEM algorithm for alignment. Requires bwa.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""


    # (FLAG) check_index: check if reference genome is already indexed
    if check_index:
        suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        if check_existence([reference + suffix for suffix in suffix_list]):
            click.echo('Index files already exist!\n Skipping reference genome indexing.')
        else:
            click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
            run(['bwa', 'index', reference])
    else:
        click.echo('Indexing reference genome %s...' % reference)
        run(['bwa', 'index', reference])


    # Align input sequences to the reference genome
    click.echo('Aligning read(s) against the reference genome...') # find way to make message fancier with custom names
    align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1, read2]
    output_name=name+'.sam'
    with open(output_name, "w+") as align_out:
        call(align_args, stdout=align_out)

    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', name + '.bam', '-T', '/tmp/lane_temp', output_name]
    run(sort_args)


    # clean up intermediary files
    if clean_up:
        click.echo('Cleaning up %s...' % output_name)
        run(['rm', name + output_name])

@align.command('bowtie2')
@click.option('--name', '-n', default='bowtie2_out', help='Name of the output file (extension will be added automatically)')
@click.option('--check_index', '-i', is_flag=True, help='Check if reference is already indexed, and if yes prompt '
                                                        'the user to skip that step.')
@click.option('--clean-up', '-c', default=False, help='Clean up intermediary files to save disk space.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(name, check_index, clean_up, reference, read1, read2):
    """Use the FM-index tool bowtie2 for alignment. Requires bowtie2.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""


    ref_basename = reference.split('.')[0] # bowtie2 uses the reference's basename (no suffix) a lot
    # (FLAG) check_index: check if reference genome is already indexed
    if check_index:
        suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        if check_existence([ref_basename + suffix for suffix in suffix_list]):
            click.echo('Index files already exist!\n Skipping reference genome indexing.')
        else:
            click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
            # Index the reference genome
            index_args = [config['bowtie2_path']+'/bowtie2-build', reference, ref_basename] # last arg is the
                                                                                            # output name
            run(index_args)
    else:
        click.echo('Indexing reference genome %s...' % reference)
        index_args = [config['bowtie2_path']+'/bowtie2-build', reference, ref_basename]
        run(index_args)


    # Align the input reads against the reference genome
    output_name = name+'.sam'
    if read2 is None: # if read is single-ended
        align_args = [config['bowtie2_path']+'/bowtie2', '-x', ref_basename, read1, '-S', output_name]
        click.echo('Aligning read %s against the reference genome...' % read1)
    else: # if paired_end
        align_args = [config['bowtie2_path']+'/bowtie2', '-x', ref_basename,
                      '-1', read1, '-2', read2, '-S', output_name]
        click.echo('Aligning reads %s %s against the reference genome...' % (read1, read2))
    run(align_args)


    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', name + '.bam', '-T', '/tmp/lane_temp', name + output_name]
    run(sort_args)


    # clean up intermediary files
    if clean_up:
        click.echo('Cleaning up %s...' % output_name)
        run(['rm', output_name])


##### POST-PROCESSING #####
@cli.command(short_help='Prepare reads for variant calling.')
@click.argument('asjfaj')
def processing(asjfaj):
    """Performs a group of steps for the post-processing of aligned reads in preparation for variant calling."""
    pass


##### VARIANT CALLING #####
@cli.group(short_help='Call variants on mapped sequences.')
def call():
    """Set of tools to call variants on files containing sequences previously aligned to a reference genome."""


@call.command('gatk')
@click.option('--name', '-n', default='gatk_out', help='Name of the output file (extension will be added automatically)')
@click.option('--platform', '-p', default='ion proton', help='Platform used in sample sequencing.')
@click.option('--library', '-l', default='library1', help='DNA preparation library identifier.')
@click.option('--sample', '-s', default='NA12878', help='The name of the sample sequenced in a read group.')
@click.option('--known-snps', '-k', default=None, type=click.Path(exists=True), help='dbSNP file containing a database \
                                                                                     of known SNPs to help improve \
                                                                                     variant calling results.')
@click.option('--clean-up', '-c', default=False, help='Clean up intermediary files to save disk space.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_gatk(name, platform, library, sample, known_snps, clean_up, reference, sample1, sample2):
    """Call variants using GATK's HaplotypeCaller.
    The GATK's HaplotypeCaller algorithm is used to call variants on aligned sequence files (samples).

    Through the aforementioned algorithm, this command calls variants on input aligned sequence files (samples),
    simplifying the operation thanks to performing all prerequisite processing steps required by GATK to call variants
    on input files.

    First, it prepares the reference by ensuring that it has been indexed through samtools faidx, and had a dictionary
    generated through Picard's CreateSequenceDictionary. Then, it prepares the input samples by creating a new file
    with proper read group (RG) information, and indexing each sample through samtools index. Lastly,
    GATK-HC is run on these prepared files.

    Specifying a dbSNP file with a list of known SNPs is highly recommended. Only one sample sequence file has to be
    specified.  Sample sequences must have been previously aligned, so that they are in the SAM/BAM format.
    A reference genome must be specified.
    """

    # samtools faidx on REFERENCE
    if check_existence([reference+'.fai']):
        click.echo('Reference faidx index file already exists!\nSkipping faidx indexing.')
    else:
        click.echo('Indexing reference file %s...' % reference)
        faidx_args = ['samtools', 'faidx', reference]
        run(faidx_args)


    # generate .dict dictionary file for REFERENCE
    dict_file = reference.split('.')[0]+'.dict'
    if check_existence([dict_file]):
        click.echo('Dictionary file %s already exists.\nSkipping reference genome dictionary file generation.' % dict_file)
    else:
        click.echo('Generating reference genome dictionary %s...' % dict_file)
        dict_vars = ['java', '-jar', config['picard_path'], 'CreateSequenceDictionary', 'R=%s' % reference,
                     'O=%s' % dict_file]
        run(dict_vars)


    # read groups
    # only works for 1 library atm
    sample_list = [sample1] + [s for s in sample2]
    sample_rg_list = [sample.split('.')[0]+'_rg.'+sample.split('.')[1] for sample in sample_list]
    if check_existence(sample_rg_list):
        click.echo('Read group information has already been added for %s.' % ', '.join(sample_list))
    else:
        for sample in sample_list:
            click.echo('Adding Read Group information for %s...' % sample)
            read_groups = {'ID': sample.split('.')[0], 'PL': platform, 'LB': library, 'PU': 'foo', 'SM': sample}
            rg_args = ['java', '-jar', config['picard_path'], 'AddOrReplaceReadGroups', 'I='+sample,
                       'O='+sample.split('.')[0]+'_rg.bam', 'RGID='+read_groups['ID'], 'RGLB='+read_groups['LB'],
                       'RGPL='+read_groups['PL'], 'RGPU='+read_groups['PU'], 'RGSM='+read_groups['SM']]
            click.echo(rg_args)
            run(rg_args)


    # samtools index on each SAMPLE
    if check_existence([sample+'.bai' for sample in sample_rg_list]):
        click.echo('Sample index .bai files already exist!\nSkipping sample indexing.')
    else:
        click.echo('Need to generate sample index .bai files!\nIndexing sample files %s...' % ', '.join(sample_rg_list))
        index_args = ['samtools', 'index'] + sample_rg_list
        run(index_args)


    # run GATK-HC
    click.echo('Calling variants on samples %s with GATK-HC...' % ', '.join(sample_rg_list))
    if known_snps is None:
        # each sample needs to be preceed by an -I so this is not working for more than one sample?
        gatk_args = ['java', '-jar', config['gatk_path'], '-R', reference, '-T', 'HaplotypeCaller', '-I'] + sample_rg_list + \
                ['-stand_call_conf', '20', '-o', name+'.vcf'] # confidence of call is not flexible atm, always 20
    else:
        gatk_args = ['java', '-jar', config['gatk_path'], '-R', reference, '-T', 'HaplotypeCaller', '-I'] + sample_rg_list + \
                    ['--dbsnp', known_snps, '-stand_call_conf', '20', '-o', name+'.vcf']
    run(gatk_args)


    # clean up intermediary files
    if clean_up:
        for sample in sample_list:
            click.echo('Cleaning up %s...' % sample)
            run(['rm', sample])