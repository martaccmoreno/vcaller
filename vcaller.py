import json
import os
from subprocess import run, getstatusoutput
import click

##### IMPORT CONFIG #####
current_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(current_dir, 'config.json'), 'r') as data_file:
    config = json.load(data_file)


##### AUXILIARY FUNCTIONS #####
def flatten_list(list_of_list):
    """Flatten a list of list into a single list."""
    return [item for sublist in list_of_list for item in sublist]


def remove_suffix(file_name):
    """Remove the suffix of a filename, e.g. 'reference.fa' becomes 'reference'"""
    return '.'.join(file_name.split('.')[:-1])


def replace_suffix(filename, new_suffix):
    """
    Replaces a filename's suffix with another user-specified suffix.
    '"""
    if new_suffix[0] == '.':
        return remove_suffix(filename) + new_suffix
    else:
        return remove_suffix(filename) + '.' + new_suffix


def check_existence(filename_list):
    """Check if files with the filenames in the list already exist in the working directory."""
    if type(filename_list) is str:
        filename_list = [filename_list]
    if sum([os.path.isfile(ifile) for ifile in filename_list]) == len(filename_list):
        return True
    else:
        return False


def tabix_index(gzipped_files):
    if type(gzipped_files) is not list: gzipped_files = [gzipped_files]
    for file in gzipped_files:
        if '.gz' in file:
            click.echo('Gunzipping %s...' % file)
            run(['gunzip', file])
            click.echo('Compressing file %s using bgzip...' % file)
            run(['bgzip', remove_suffix(file)])
            click.echo('Indexing file %s using tabix....' % file)
            run(['tabix', file])
        else:
            click.echo('Compressing file %s using bgzip...' % file)
            run(['bgzip', file])
            click.echo('Indexing file %s using tabix....' % file)
            run(['tabix', file + '.gz'])


def cleanup(files_or_dirs):
    """
    Permanently remove one or more files or directories.
    """
    if type(files_or_dirs) is not list: files_or_dirs = [files_or_dirs]
    click.echo('Removing the following files/directories: %s' % ', '.join(files_or_dirs))
    for file_or_dir in files_or_dirs:
        if os.path.isfile(file_or_dir):
            run(['rm', file_or_dir])
        elif os.path.isdir(file_or_dir):
            run(['rm', '-r', file_or_dir])
        else:
            click.echo('Invalid input: %s is neither a file nor a directory...' % file_or_dir)


def bed_intersect(vcf, bed, out=None, clean=False):
    """
    Intersect a vcf file with a bed file, obtaining a second vcf file with only the regions defined in the bed file.
    """
    if out is None:
        out = remove_suffix(vcf)+'_exome.vcf'
    with open(out, 'w+') as out:
        click.echo('Intersecting vcf %s with bed file regions in %s...\n' % (vcf, bed))
        run('bedtools', 'intersect', '-header', '-a', vcf, '-b', bed, stdout=out)
    if clean:
        cleanup(vcf)


### BROADCAST FUNCTIONS
def broadcast_ref_index(suffixes, reference):
    if check_existence(suffixes):
        click.echo('The following index files already exist:\n %s' % ' '.join(suffixes))
        click.echo('Skipping reference genome indexing...\n')
    else:
        click.echo('Need to generate index files for %s!\n Indexing reference genome %s...\n' % reference)

def broadcast_alignment(reads, reference):
    if len(reads) == 1:
        click.echo('Aligning read %s against the reference genome %s...\n' % (reads[0]. reference))
    else:
        click.echo('Aligning read(s) %s against the reference genome %s...\n' % (' '.join(reads), reference))


##### MAIN GROUP #####
@click.group()
@click.version_option()
def cli():
    """
    Vcaller, a CLI capable of evoking multiple pre-existing bioinformatics tools, ecletically grouping them into
    commands and subcommands that will perform common variant calling and benchmarking routines.
    """


##### SEQUENCE ALIGNMENT #####
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
def align_bowtie2(output, no_clean, reference, read1, read2):
    """Use the FM-index tool Bowtie 2 for alignment.
    Requires Bowtie 2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    suffixes = [remove_suffix(reference) + suffix for suffix in suffix_list]
    broadcast_ref_index(suffixes, reference)
    if not check_existence(suffixes):
        index_args = [config['filePaths']['bowtie2'] + '/bowtie2-build', reference, remove_suffix(reference)]
        run(index_args)

    sam_output = replace_suffix(output, 'sam')
    broadcast_alignment([read1, read2], reference)
    align_args = [config['filePaths']['bowtie2'] + '/bowtie2', '-x', remove_suffix(reference), '-S', sam_output,
                  read1]
    if read2 is not None:
        align_args += [read2]
    run(align_args)

    click.echo('Sorting and converting %s to the BAM format...\n' % sam_output)
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', os.path.join('/tmp/', replace_suffix(
        os.path.basename(output)), 'tmp'), sam_output, sam_output]
    run(sort_args)

    if no_clean is False:
        cleanup(sam_output)


@align.command('bwa')
@click.option('--output', '-o', default='bwa_out.bam', help='Name of the output file.')
@click.option('--nthreads', '-t', default='1', help='Number of CPU threads to use during the alignment step.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('read1', type=click.Path(exists=True))
@click.argument('read2', required=False, type=click.Path(exists=True))
def align_bwa(output, nthreads, no_clean, reference, read1, read2):
    """Use the BWA-MEM algorithm for alignment. Requires bwa.
    It is only mandatory to include the reference genome file and a sample read as arguments.
    If dealing with paired-end reads, a second sequence file containing the second mate-pair read may be included."""

    # Check if reference genome is already indexed
    suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    if check_existence([reference + suffix for suffix in suffix_list]):
        click.echo('Index files already exist!\n Skipping reference genome indexing.\n')
    else:
        click.echo('Need to generate index files!\n Indexing reference genome %s...' % reference)
        run(['bwa', 'index', reference])

    # Align input sequences to the reference genome
    sam_output = replace_suffix(output, 'sam')
    if check_existence([sam_output]):
        click.echo('Aligned SAM read file already exists!')
    else:
        click.echo('Aligning reads against the reference genome...')
        align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1]
        if read2 is not None:
            align_args += [read2]
        with open(sam_output, "w+") as align_out:
            run(align_args, stdout=align_out)

    click.echo('Sorting and converting to BAM...')
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', os.path.join('/tmp/', replace_suffix(
        os.path.basename(output)), 'tmp'), sam_output]
    run(sort_args)

    # Remove intermediary files
    click.echo('Cleaning up %s...' % sam_output)
    run(['rm', sam_output])


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

    suffix_list = ['.tmap.anno', '.tmap.bwt', '.tmap.pac', '.tmap.sa']
    suffixes = [reference + suffix for suffix in suffix_list]
    if check_existence(suffixes):
        click.echo('Index files already exist! Skipping reference genome indexing.')
    else:
        click.echo('Need to generate index files! Indexing reference genome %s...' % reference)
        index_args = [config['filePaths']['tmap'], 'index', '-f', reference]
        run(index_args)

    if read2 is None:  # if read is single-ended
        align_args = [config['filePaths']['tmap'], 'map1', '-o', '2', '-f', reference, '-r', read1]
        if read2 is not None:
            align_args += [read2]
        if 'gz' in read1.split('.') or 'gz' in read2.split('.'):
            align_args += ['--input-gz']
        click.echo('Aligning reads against the reference genome...')
        with open(output, "w+") as align_out:
            run(align_args, stdout=align_out)


##### VARIANT CALLING #####
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
def call_bcftools(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2):
    """Call variants using SAMtools's BCFtools.

    This command calls variants on input aligned sequence files (samples) after calculating their genotype likelihoods.

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """

    sample_list = [sample1] + [s for s in sample2]
    bcf_output = replace_suffix(output, 'bcf')

    click.echo('Calculating genotype likelihoods for %s...' % ', '.join(sample_list))
    if count_orphans:
        mpileup_args = ['bcftools', 'mpileup', '-AOb', '-o', bcf_output, '-f', reference] + sample_list
    else:
        mpileup_args = ['bcftools', 'mpileup', '-Ob', '-o', bcf_output, '-f', reference] + sample_list
    run(mpileup_args)

    click.echo('Calling variants on %s with BCFtools...' % ', '.join(sample_list))
    call_args = ['bcftools', 'call', '-vmO', 'v', '-o', output, bcf_output]
    run(call_args)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)

    if no_clean is False:
        cleanup(bcf_output)


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

    sample_list = [sample1] + [s for s in sample2]

    click.echo('Calling variants on %s using Freebayes...' % ', '.join(sample_list))
    call_args = [config['filePaths']['freebayes'], '-f', reference] + sample_list
    with open(output, 'w+') as call_out:
        run(call_args, stdout=call_out)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)


@call.command('gatk')
@click.option('--output', '-o', default='gatk_out.vcf',
              help='Name of the output file (extension will be added automatically)')
@click.option('--dbsnp', default=None, type=click.Path(exists=True), help='dbSNP file containing a database of '
                                                                          'known SNP IDs.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_gatk(output, dbsnp, exome_regions, reference, sample1, sample2):
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

    if check_existence([reference + '.fai']):
        click.echo('Reference faidx index file already exists! Skipping faidx indexing.')
    else:
        click.echo('Indexing reference file %s...' % reference)
        run(['samtools', 'faidx', reference])

    dict_file = replace_suffix(reference, 'dict')
    if check_existence([dict_file]):
        click.echo(
            'Dictionary file %s already exists. Skipping reference genome dictionary file generation.' % dict_file)
    else:
        click.echo('Generating reference genome dictionary %s...' % dict_file)
        dict_vars = ['java', '-jar', config['filePaths']['picard'], 'CreateSequenceDictionary', 'R=%s' % reference,
                     'O=%s' % dict_file]
        run(dict_vars)

    sample_list = [sample1] + [s for s in sample2]
    for smpl in sample_list:
        if check_existence([replace_suffix(smpl, 'bai')]) or check_existence(smpl + '.bai'):
            click.echo('Sample index .bai files already exist! Skipping sample indexing.')
        else:
            click.echo('Need to generate sample index .bai file!\nIndexing sample file %s...' % smpl)
            run(['samtools', 'index', smpl])

    click.echo('Calling variants on samples %s with GATK-HC...' % ', '.join(sample_list))
    if dbsnp is None:
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference] + \
                    flatten_list([['-I'] + [sample_list[i]] for i in range(len(sample_list))]) + ['-O', output]
    else:
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference] + \
                    flatten_list([['-I'] + [sample_list[i]] for i in range(len(sample_list))]) + \
                    ['--dbsnp', dbsnp, '-O', output]
    run(gatk_args)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)


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

    click.echo('Calling variants on %s with TVC...' % sample)
    call_args = [config['filePaths']['tvc'], '-i', sample, '-r', reference, '-o', output_dir]
    if exome_regions:
        call_args += ['-b', exome_regions]
    run(call_args)


@call.command('varscan2')
@click.option('--output', '-o', default='bcftools_out.vcf', help='Name of the output file.')
@click.option('--count-orphans', '-A', is_flag=True, help='Count reads with anomolous mate pairs.')
@click.option('--exome-regions', '-e', default=None, help='Bed file to restrict output regions to the exome.')
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('reference', type=click.Path(exists=True))
@click.argument('sample1', type=click.Path(exists=True))
@click.argument('sample2', required=False, type=click.Path(exists=True), nargs=-1)
def call_tvc(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2):
    """Call variants using Varscan2.

    Only one sample sequence file has to be specified. Sample sequences must have been previously aligned, so that they
    are in the SAM/BAM format. A reference genome must be provided.
    """

    sample_list = [sample1] + [s for s in sample2]

    click.echo('Creating mpipleup file using %s...' % ', '.join(sample_list))
    mpileup_file = replace_suffix(output, 'pileup')
    if count_orphans:
        pileup_args = ['samtools', 'mpileup', '-A', '-f', reference] + sample_list
    else:
        pileup_args = ['samtools', 'mpileup', '-f', reference] + sample_list
    with open(mpileup_file, 'w+') as pileup_out:
        run(pileup_args, stdout=pileup_out)

    click.echo('Calling variants on %s with Varscan2...' % mpileup_file)
    call_args = ['java', '-jar', config['filePaths']['varscan2'], 'mpileup2cns', mpileup_file, '--output-vcf', '1',
                 '--variants', '1', '--p-value', '0.10', '--min-coverage', '2']
    with open(output, 'w+') as call_out:
        run(call_args, stdout=call_out)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)

    if no_clean is False:
        cleanup(mpileup_file)


##### POST-PROCESSING #####
@cli.command('process', short_help='Prepare reads for variant calling.')
@click.option('--output-name', '-o', default=None,
              help='Name of the output file (extension will be added automatically)')
@click.option('--output-dir', '-d', default=None,
              help='Name of output directory; by default save in the same directory as the final output.')
@click.option('--readgroup-info', default=None, type=str, help='Add read group information  to the sample, which MUST '
                                                               'follow the format below:\n'
                                                               'ID:identifier,PU:platform_unit,' '\n'
                                                               'PL:platform,SM:sample,LB:library' '\n')
@click.option('--add-known-snps', '-s', default='', help='Additional files containing known SNP information.',
              multiple=True)
@click.option('--add-known-indels', '-i', default='', help='Additional files containing known indel information.',
              multiple=True)
@click.option('--no-clean', is_flag=True, help='Do not remove intermidiary files')
@click.argument('known-indels', required=True,
                type=click.Path(exists=True))  # ADD OPTIONS FOR MORE KNOWN INDELS AND SNPS
@click.argument('known-snps', required=True, type=click.Path(exists=True))
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('sample', required=True, type=click.Path(exists=True))
def process(output_name, output_dir, readgroup_info, add_known_snps, add_known_indels, no_clean, known_indels, known_snps,
            reference,
            sample):
    """Performs a group of steps for the post-processing in preparation for variant calling
    on one SAM/BAM sampl file. A must do for running the gatk subcommand under call."""

    if output_dir is None:
        output_dir = os.path.dirname(output_name)
    smpl_name, smpl_extension = remove_suffix(os.path.basename(sample)), '.' + sample.split('.')[-1]

    if not check_existence(reference + '.fai'):
        click.echo('Generating faidx index for %s...' % reference)
        run(['samtools', 'faidx', reference])
    if not check_existence([remove_suffix(reference) + '.dict']):
        click.echo('Generating sequence dictionary for %s...' % reference)
        run([config['filePaths']['gatk4'], 'CreateSequenceDictionary', '-R', reference])

    if (smpl_extension is not '.bam') or (getstatusoutput('samtools index ' + sample)[0] != 0):
        sorted_output = replace_suffix(sample, 'bam')
        if not check_existence([sorted_output]):
            click.echo('Sorting and converting %s to BAM...' % sample)
            sort_args = ['samtools', 'sort', '-O', 'bam', '-o', sorted_output, '-T',
                         os.path.join('/tmp/', smpl_name + '.temp'), sample]
            run(sort_args)
            click.echo('Indexing %s...' % sorted_output)
            run(['samtools', 'index', sorted_output])
        sample = sorted_output

    if check_existence(sample + '.bai'):  # to avoid unnecessary 'file not found' errors
        run(['rm', sample + '.bai'])

    # dedupping
    dup_output = os.path.join(output_dir, smpl_name + '.DUP.bam')
    dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', sample,
                '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                '-M', os.path.join(output_dir, smpl_name + '.metrics')]
    if readgroup_info:
        # More info on RGs: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        rg_output = os.path.join(output_dir, smpl_name + '.RG.bam')
        if not check_existence(rg_output):
            click.echo('Adding Read Group information to %s...' % sample)
            rg_info = readgroup_info.split(',')
            read_groups = {'ID': rg_info[0].split(':')[1], 'PU': rg_info[1].split(':')[1],
                           'PL': rg_info[2].split(':')[1], 'SM': rg_info[3].split(':')[1],
                           'LB': rg_info[4].split(':')[1]}
            rg_args = [config['filePaths']['gatk4'], 'AddOrReplaceReadGroups', '-I', sample,
                       '-O', rg_output, '-RGID', read_groups['ID'], '-RGLB', read_groups['LB'],
                       '-RGPL', read_groups['PL'].upper(), '-RGPU', read_groups['PU'], '-RGSM', read_groups['SM'],
                       '-SO', 'coordinate']
            run(rg_args)
        dup_output = replace_suffix(rg_output, '.DUP') + '.bam'
        dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', rg_output,
                    '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                    '-M', os.path.join(output_dir, smpl_name + '.metrics')]

    if not check_existence(dup_output):
        click.echo('Marking and removing duplicates for %s...' % sample)
        run(dup_args)

    if not check_existence(dup_output+'.bai'):
        run(['samtools', 'index', dup_output])

    # Realign around indels; using gatk3 because of this step
    # https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
    intervals_output = os.path.join(output_dir, smpl_name + '.intervals')
    if not check_existence([intervals_output]):
        click.echo('Creating indel realignment intervals for %s...' % sample)
        intervals_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'RealignerTargetCreator', '-R', reference,
                          '-I', dup_output, '-o', intervals_output, '--known',
                          known_indels]
        run(intervals_args)
    realign_output = replace_suffix(dup_output, 'RLGN') + '.bam'
    if not check_existence([realign_output]):
        click.echo('Applying indel realignment based on the intervals for %s...' % sample)
        realign_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'IndelRealigner', '-R', reference,
                        '-I', dup_output, '-targetIntervals', intervals_output, '-known', known_indels] + \
                       flatten_list([['-known'] + [add_known_indels[i]] for i in range(len(add_known_indels))]) + \
                       ['-o', realign_output]
        run(realign_args)

    # BQSR
    table_output = os.path.join(output_dir, smpl_name + '.table')
    if not check_existence([table_output]):
        click.echo('Creating base score recalibration table for %s...' % sample)
        table_args = [config['filePaths']['gatk4'], 'BaseRecalibrator', '-R', reference,
                      '--known-sites', known_snps] + \
                     flatten_list([['--known-sites'] + [add_known_snps[i]] for i in range(len(add_known_snps))]) + \
                     ['-I', realign_output, '-O', table_output]
        run(table_args)
    if output_name is None:
        bqsr_output = os.path.join(output_dir, smpl_name + '.processed.bam')
    else:
        bqsr_output = output_name
    if not check_existence([bqsr_output]):
        click.echo('Running base score recalibration on %s...' % sample)
        bqsr_args = [config['filePaths']['gatk4'], 'ApplyBQSR',
                     '-I', realign_output, '-bqsr', table_output, '-O', bqsr_output]
        run(bqsr_args)

    # clean up intermediary files -- but only after we have the final file
    # still not working 100% -- IT STILL RUNS EVEN IF BSQR WASN'T SUCCESSFUL
    if check_existence([bqsr_output]) and no_clean is False:
        files_to_rmv = [dup_output, intervals_output, realign_output, table_output,
                    os.path.join(output_dir, smpl_name + '.metrics')]
        if readgroup_info:
            files_to_rmv += [rg_output]
        index_files = [replace_suffix(item, 'bai') for item in [file for file in files_to_rmv if '.bam' in file]]
        cleanup(files_to_rmv+index_files)


##### COMPARE #####
@cli.command('compare', short_help='Compare two sets of called variants, '
                                   'with one of them being assumed to be the baseline set.')
@click.option('--output-dir', '-o', default='comparison', # change to just output
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
def compare(output_dir, bed_file, evaluation_regions, score_field, sample, no_clean, reference, baseline, calls):
    """
    Evokes rtgtool's vcfeval to compare a set of baseline calls against a set of query calls, outputting a GA4GH-compliant
    annotated VCF. Next, this VCF is passed into hap.py's qfy.py method in order to compute metrics, namely raw counts
    of TP/FP/FN, as well as their associated precision and recall.
    Because of the way vcfeval's algorithm works, a SDF format file of the reference genome will have to be generated.
    Optionally, regions of the genome wherein to produce the comparison (e.g. exome capture kit regions), as well as
    those considered to be of high confidence, can be defined by providing a BED file.
    To select a pair of samples from each variant data set, use the format <baseline_sample>,<calls_sample>.
    """

    # The reference genome must be converted to SDF
    sdf_ref = os.path.join(replace_suffix(reference, 'sdf'))
    if os.path.isdir(sdf_ref) is False:
        click.echo('Converting reference genome %s to the SDF format...' % reference)
        fastq2sdf_args = [config['filePaths']['rtg'], 'format', '-o', sdf_ref, reference]
        run(fastq2sdf_args)
    else:
        click.echo('Reference genome %s has already been convert to the SDF format under %s' %
                   (os.path.basename(reference), sdf_ref))

    # Check if baseline and calls are tabix-indexed as this is a must for using vcfeval
    if not check_existence(baseline + '.tbi'):
        tabix_index(baseline)
        baseline += '.gz'
    if not check_existence(calls + '.tbi'):
        tabix_index(calls)
        calls += '.gz'

    click.echo('Creating directory %s...' % output_dir)

    # Create GA4GH-compliant annotated VCFs
    rtg_out = os.path.join(output_dir, os.path.split(output_dir)[-1] + '-vcfeval')
    if not os.path.isdir(rtg_out):
        click.echo('Comparing baseline %s against call set %s using vcfeval...' % (baseline, calls))
        rtg_args = [config['filePaths']['rtg'], 'vcfeval', '-o', rtg_out, '--vcf-score-field',
                    score_field, '--template', sdf_ref, '--baseline', baseline, '--calls', calls, '-m', 'ga4gh']
        if bed_file is not None: rtg_args += ['--bed-regions', bed_file]
        if evaluation_regions is not None: rtg_args += ['--evaluation-regions', evaluation_regions]
        if sample is not None: rtg_args += ['--sample', sample]
        run(rtg_args)

    click.echo('Moving to directory %s...' % output_dir)
    initial_path = os.getcwd()
    os.chdir(output_dir)
    click.echo('Running qfy.py on %s...' % os.path.join(rtg_out, 'output.vcf.gz'))
    qfy_args = [config['filePaths']['qfy.py'], '-t', 'ga4gh', '--verbose', '--adjust-conf-regions',
                os.path.normpath(os.path.join(initial_path, evaluation_regions)), '--reference',
                os.path.normpath(os.path.join(initial_path, reference)), '-o', os.path.split(output_dir)[-1],
                '--write-vcf', '--write-counts',
                os.path.join(os.path.basename(rtg_out), 'output.vcf.gz')]
    run(qfy_args)
    click.echo('Returning to %s...' % initial_path)
    os.chdir(initial_path)

    if no_clean is False:
        cleanup(rtg_out)