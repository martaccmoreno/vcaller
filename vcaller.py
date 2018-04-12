import json
import os
from subprocess import run

import click

##### IMPORT CONFIG #####
# find where the script directory is (=/= working directory)
current_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(current_dir, 'config.json')) as data_file:
    config = json.load(data_file)


##### AUXILIARY FUNCTIONS #####
def check_existence(filename_list):
    """Check if files with the filenames in the list already exist in the working directory."""
    if sum([os.path.isfile(ifile) for ifile in filename_list]) == len(filename_list):
        return True
    else:
        return False


def flatten_list(list_of_list):
    "Flatten a list of list into a single list."
    return [item for sublist in list_of_list for item in sublist]


##### MAIN GROUP #####
@click.group()
@click.version_option()
def cli():
    """Vcaller
    One-stop application that allows the user to leverage several existing bioinformatics tools that constitute variant
    calling pipelines. If desired, the recommended settings for each command can be customized.
    """


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
    output_name = ''.join(output.split('.')[:-1])
    sam_output = output_name + '.sam'
    if check_existence([sam_output]):
        click.echo('Aligned SAM read file already exists!')
    else:
        click.echo(
            'Aligning read(s) against the reference genome...')  # find way to make message fancier with custom names
        align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1]
        if read2 is not None:
            align_args += [read2]
        with open(sam_output, "w") as align_out:
            run(align_args, stdout=align_out)

    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', '/tmp/' + output_name + '_temp', sam_output]
    run(sort_args)

    # clean up intermediary files
    click.echo('Cleaning up %s...' % sam_output)
    run(['rm', sam_output])


@align.command('bowtie2')
@click.option('--output', '-o', default='bowtie2_output.bam',
              help='Name of the output file (extension will be added automatically)')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(output, reference, read1, read2):
    """Use the FM-index tool bowtie2 for alignment. Requires bowtie2.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    if check_existence([''.join(reference.split('.')[:-1]) + suffix for suffix in suffix_list]):
        click.echo('Index files already exist!\n Skipping reference genome indexing.')
    else:
        click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
        index_args = [config['filePaths']['bowtie2'] + '/bowtie2-build', reference, ''.join(reference.split('.')[:-1])]
        run(index_args)

    output_basename = os.path.basename(''.join(output.split('.')[:-1]))
    sam_output = output_basename + '.sam'
    if read2 is None:  # if read is single-ended
        align_args = [config['filePaths']['bowtie2'] + '/bowtie2', '-x', ''.join(reference.split('.')[:-1]), read1,
                      '-S', sam_output]
        click.echo('Aligning read %s against the reference genome...' % read1)
    else:  # if paired_end
        align_args = [config['filePaths']['bowtie2'] + '/bowtie2', '-x', ''.join(reference.split('.')[:-1]),
                      '-1', read1, '-2', read2, '-S', sam_output]
        click.echo('Aligning reads %s %s against the reference genome...' % (read1, read2))
    run(align_args)

    # Sort and convert to BAM
    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', '/tmp/' + output_basename + '_temp', sam_output]
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
def call_gatk(output, known_snps, reference, sample1, sample2):
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
        dict_vars = ['java', '-jar', config['filePaths']['picard'], 'CreateSequenceDictionary', 'R=%s' % reference,
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
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference, '-I'] + sample_list + \
                    ['-O', output]
    else:
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference, '-I'] + sample_list + \
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
@click.option('--output-dir', '-d', default='.',
              help='Name of output directory; by default save to current directory.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_tvc(output_dir, reference, sample1, sample2):
    """Call variants using TorrentVariantCaller (TVC).

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """

    sample_list = [sample1] + [s for s in sample2]
    click.echo('Calling variants on %s with TVC...' % ', '.join(sample_list))
    call_args = [config['filePaths']['tvc'], '-i', ','.join(sample_list), '-r', reference, '-o', output_dir]
    run(call_args)


##### POST-PROCESSING #####
@cli.command('process', short_help='Prepare reads for variant calling.')
@click.option('--output-name', '-o', default=None,
              help='Name of the output file (extension will be added automatically)')
@click.option('--output-dir', '-d', default='',
              help='Name of output directory; by default save to current directory.')
@click.option('--readgroup-info', default=None, type=str, help='Add read group information  to the sample, which MUST '
                                                               'follow the format below:\n'
                                                               r'\tID:identifier\tPU:platform_unit' '\n'
                                                               r'\tPL:platform\tSM:sample\tLB:library' '\n')
@click.option('--add-known-snps', '-s', default='', help='Additional files containing known SNP information.',
              multiple=True)
@click.option('--add-known-indels', '-i', default='', help='Additional files containing known indel information.',
              multiple=True)
@click.argument('known-indels', required=True,
                type=click.Path(exists=True))  # ADD OPTIONS FOR MORE KNOWN INDELS AND SNPS
@click.argument('known-snps', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('sample', required=True, type=click.Path(exists=True))
def process(output, output_dir, readgroup_info, add_known_snps, add_known_indels, known_indels, known_snps, reference,
            sample):
    """Performs a group of steps for the post-processing in preparation for variant calling
    on one SAM/BAM sampl file. A must do for running the gatk subcommand under call."""

    if output is None:
        smpl_name, smpl_extension = '.'.join(sample.split('.')[:-1]), sample.split('.')[-1]
        smpl_name = os.path.basename(smpl_name)
        smpl_extension = '.' + smpl_extension
    else:
        smpl_name = os.path.basename(output)
        smpl_extension = output.split('.')[-1]

    # sort and convert SAM extension files to BAM
    if smpl_extension.lower() is '.sam':
        click.echo('Sorting and converting %s to BAM...' % sample)
        sort_args = ['samtools', 'sort', '-O', 'bam', '-o', smpl_name + '.bam', '-T', '/tmp/lane_temp', sample]
        run(sort_args)
        smpl_extension = '.bam'

    # More info on RGs: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    if readgroup_info is not None:
        rg_output = output_dir + smpl_name + '.RG' + smpl_extension
        if not check_existence([rg_output]):
            click.echo('Adding Read Group information to %s...' % sample)
            rg_info = readgroup_info.split(r'\t')
            click.echo(rg_info)
            read_groups = {'ID': rg_info[0].split(':')[1], 'PU': rg_info[1].split(':')[1],
                           'PL': rg_info[2].split(':')[1], 'SM': rg_info[3].split(':')[1],
                           'LB': rg_info[4].split(':')[1]}
            click.echo(read_groups)
            rg_args = [config['filePaths']['gatk4'], 'AddOrReplaceReadGroups', '-I', sample,
                       '-O', rg_output, '-RGID', read_groups['ID'], '-RGLB', read_groups['LB'],
                       '-RGPL', read_groups['PL'].upper(), '-RGPU', read_groups['PU'], '-RGSM', read_groups['SM']]
            run(rg_args)
        click.echo('Marking and removing duplicates for %s...' % smpl_name)
        dup_output = '.'.join(rg_output.split('.')[:-1]) + '.DUP' + smpl_extension
        dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', rg_output,
                    '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                    '-M', smpl_name + '.metrics']
    else:
        click.echo('Marking and removing duplicates for %s...' % smpl_name)
        dup_output = output_dir + smpl_name + '.DUP' + smpl_extension
        dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', sample,
                    '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                    '-M', smpl_name + '.metrics']
    if not check_existence([dup_output]):
        run(dup_args)

    # Realign around indels
    # Using gatk3 because of this
    # https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
    intervals_output = smpl_name + '.intervals'
    if not check_existence([intervals_output]):
        click.echo('Creating indel realignment intervals for %s...' % smpl_name)
        intervals_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'RealignerTargetCreator', '-R', reference,
                          '-I', dup_output, '-o', intervals_output, '--known',
                          known_indels]  # only one set of known atm
        run(intervals_args)
    realign_output = '.'.join(dup_output.split('.')[:-1]) + '.RLGN' + smpl_extension
    if not check_existence([realign_output]):
        click.echo('Applying indel realignment based on the intervals for %s...' % smpl_name)
        realign_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'IndelRealigner', '-R', reference,
                        '-I', dup_output, '-targetIntervals', intervals_output, '-known', known_indels] + \
                       flatten_list([['-known'] + [add_known_indels[i]] for i in range(len(add_known_indels))]) + \
                       ['-o', realign_output]
        run(realign_args)

    # BQSR
    table_output = smpl_name + '.table'  # should be the ACTUAL name for the file...
    if not check_existence([table_output]):
        click.echo('Creating base score recalibration table for %s...' % smpl_name)
        table_args = [config['filePaths']['gatk4'], 'BaseRecalibrator', '-R', reference,
                      '--known-sites', known_snps] + \
                     flatten_list([['--known-sites'] + [add_known_snps[i]] for i in range(len(add_known_snps))]) + \
                     ['-I', realign_output, '-O', table_output]
        run(table_args)
    bqsr_output = output_dir + smpl_name + '.processed' + smpl_extension
    if not check_existence([bqsr_output]):
        click.echo('Running base score recalibration on %s...' % smpl_name)
        bqsr_args = [config['filePaths']['gatk4'], 'ApplyBQSR',
                     '-I', realign_output, '-bqsr', table_output, '-O', bqsr_output]
        run(bqsr_args)

    # clean up intermediary files -- but only after we have the final file
    if check_existence([bqsr_output]) and readgroup_info is not None:
        click.echo('Cleaning up %s...' % ', '.join([rg_output, dup_output, intervals_output, realign_output,
                                                    table_output, smpl_name + '.metrics']))
        run(['rm', rg_output, dup_output, intervals_output, realign_output, table_output])
    elif check_existence([bqsr_output]):
        click.echo('Cleaning up %s...' % ', '.join([dup_output, intervals_output, realign_output,
                                                    table_output, smpl_name + '.metrics']))
    run(['rm', dup_output, intervals_output, realign_output, table_output])


##### PREP #####

@cli.command('tabix', short_help='Quickly tabix gunzipped known files.')
@click.argument('gzipped_files', required=True, type=click.Path(exists=True), nargs=-1)
def tabix(gzipped_files):
    for file in gzipped_files:
        click.echo('Gunzipping %s...' % file)
        run(['gunzip', file])
        click.echo('Compressed gunzipped file %s using bgzip...' % file)
        run(['bgzip', ''.join(file.split('.')[:-1])])
        click.echo('Indexing file %s using tabix....' % file)
        run(['tabix', file])


##### FILTERING #####
@cli.command('filter', short_help='Filter variants in a vcf file.')
@click.argument('reference', required=True)
@click.argument('variants', required=True, type=click.Path(exists=True), nargs=-1)
def process(reference, variants):
    """Desc"""
    variant_list = list(variants)
