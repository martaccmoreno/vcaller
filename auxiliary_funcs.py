import click
import os
import subprocess
import time


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
            click.echo('\nGunzipping %s...' % file)
            subprocess.run(['gunzip', file])
            click.echo('\nCompressing file %s using bgzip...' % file)
            subprocess.run(['bgzip', remove_suffix(file)])
            click.echo('\nIndexing file %s using tabix...' % file)
            subprocess.run(['tabix', file])
        else:
            click.echo('\nCompressing file %s using bgzip...' % file)
            subprocess.run(['bgzip', file])
            click.echo('\nIndexing file %s using tabix....' % file)
            subprocess.run(['tabix', file + '.gz'])


def cleanup(files_or_dirs):
    """
    Permanently remove one or more files or directories.
    """
    if type(files_or_dirs) is not list: files_or_dirs = [files_or_dirs]
    click.echo('\nRemoving the following files/directories: %s' % ', '.join(files_or_dirs))
    for file_or_dir in files_or_dirs:
        if os.path.isfile(file_or_dir):
            subprocess.run(['rm', file_or_dir])
        elif os.path.isdir(file_or_dir):
            subprocess.run(['rm', '-r', file_or_dir])
        else:
            click.echo('\nInvalid input: %s is neither a file nor a directory...' % file_or_dir)


def bed_intersect(vcf, bed, out=None, clean=False):
    """
    Intersect a vcf file with a bed file, obtaining a second vcf file with only the regions defined in the bed file.
    """
    if out is None:
        out = remove_suffix(vcf)+'_exome.vcf'
    with open(out, 'w+') as out:
        click.echo('\nIntersecting vcf %s with bed file regions in %s...' % (vcf, bed))
        subprocess.run('bedtools', 'intersect', '-header', '-a', vcf, '-b', bed, stdout=out)
    if clean:
        cleanup(vcf)

def timestamp():
    return time.strftime("[%Y-%m-%d %H:%M]")

###########################
### BROADCAST FUNCTIONS ###
###########################

def broadcast_ref_index(suffixes, reference):
    if check_existence(suffixes):
        click.echo('\n %s The following index files already exist:\n %s' % (timestamp(), ' '.join(suffixes)))
        click.echo('\n %s Skipping reference genome indexing...' % timestamp())
    else:
        click.echo('\n %s Need to generate index files for %s!'
                   '\n %s Indexing reference genome %s...' % (timestamp(), reference, timestamp(), reference))


def broadcast_alignment(reads, reference):
    if len(reads) == 1:
        click.echo('\n %s Aligning read %s against the reference genome %s...' % (timestamp(), reads[0], reference))
    else:
        click.echo('\n %s Aligning read(s) %s against the reference genome %s...' % (timestamp(), ' '.join(reads),
                                                                                     reference))