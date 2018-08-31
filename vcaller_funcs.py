import click
from subprocess import run
import os
import json


#####################
### IMPORT CONFIG ###
#####################

current_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(current_dir, 'config.json'), 'r') as data_file:
    config = json.load(data_file)


###########################
### AUXILIARY FUNCTIONS ###
###########################

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
            click.echo('\nGunzipping %s...\n' % file)
            run(['gunzip', file])
            click.echo('\nCompressing file %s using bgzip...\n' % file)
            run(['bgzip', remove_suffix(file)])
            click.echo('\nIndexing file %s using tabix....\n' % file)
            run(['tabix', file])
        else:
            click.echo('\nCompressing file %s using bgzip...\n' % file)
            run(['bgzip', file])
            click.echo('\nIndexing file %s using tabix....\n' % file)
            run(['tabix', file + '.gz'])


def cleanup(files_or_dirs):
    """
    Permanently remove one or more files or directories.
    """
    if type(files_or_dirs) is not list: files_or_dirs = [files_or_dirs]
    click.echo('\nRemoving the following files/directories: %s\n' % ', '.join(files_or_dirs))
    for file_or_dir in files_or_dirs:
        if os.path.isfile(file_or_dir):
            run(['rm', file_or_dir])
        elif os.path.isdir(file_or_dir):
            run(['rm', '-r', file_or_dir])
        else:
            click.echo('\nInvalid input: %s is neither a file nor a directory...\n' % file_or_dir)


def bed_intersect(vcf, bed, out=None, clean=False):
    """
    Intersect a vcf file with a bed file, obtaining a second vcf file with only the regions defined in the bed file.
    """
    if out is None:
        out = remove_suffix(vcf)+'_exome.vcf'
    with open(out, 'w+') as out:
        click.echo('\nIntersecting vcf %s with bed file regions in %s...\n' % (vcf, bed))
        run('bedtools', 'intersect', '-header', '-a', vcf, '-b', bed, stdout=out)
    if clean:
        cleanup(vcf)


### BROADCAST FUNCTIONS ###
def broadcast_ref_index(suffixes, reference):
    if check_existence(suffixes):
        click.echo('\nThe following index files already exist:\n %s' % ' '.join(suffixes))
        click.echo('\nSkipping reference genome indexing...\n')
    else:
        click.echo('\nNeed to generate index files for %s!\nIndexing reference genome %s...\n' % (reference, reference))


def broadcast_alignment(reads, reference):
    if len(reads) == 1:
        click.echo('\nAligning read %s against the reference genome %s...\n' % (reads[0], reference))
    else:
        click.echo('\nAligning read(s) %s against the reference genome %s...\n' % (' '.join(reads), reference))


######################
### MAIN FUNCTIONS ###
######################

### ALIGNERS ####
def func_align_bowtie2(output, no_clean, reference, read1, read2=''):
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

    click.echo('\nSorting and converting %s to the BAM format...\n' % sam_output)
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T', os.path.join('/tmp/', replace_suffix(
        os.path.basename(output)), 'tmp'), sam_output, sam_output]
    run(sort_args)

    if no_clean is False:
        cleanup(sam_output)

def func_align_bwa(output, nthreads, no_clean, reference, read1, read2=''):
    suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    suffixes = [remove_suffix(reference) + suffix for suffix in suffix_list]
    broadcast_ref_index(suffixes, reference)
    if not check_existence([reference + suffix for suffix in suffix_list]):
        run(['bwa', 'index', reference])

    sam_output = replace_suffix(output, 'sam')
    if check_existence([sam_output]):
        click.echo('\nAligned reads SAM file already exists!\n')
    else:
        broadcast_alignment([read1, read2], reference)
        align_args = ['bwa', 'mem', '-M', '-t', nthreads, reference, read1]
        if read2 is not None:
            align_args += [read2]
        with open(sam_output, "w+") as align_out:
            run(align_args, stdout=align_out)

    click.echo('\nSorting and converting %s to the BAM format...\n' % sam_output)
    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                 os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
    run(sort_args)

    if no_clean is False:
        cleanup(sam_output)

def func_align_tmap(output, reference, read1, read2=''):
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