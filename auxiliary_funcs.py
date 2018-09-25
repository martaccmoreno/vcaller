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
        subprocess.run(['bedtools', 'intersect', '-header', '-a', vcf, '-b', bed], stdout=out)
    if clean:
        cleanup(vcf)

def timestamp():
    """
    Prints a pretty timestamp.
    """
    return time.strftime("[%Y-%m-%d %H:%M]")

###########################
### BROADCAST FUNCTIONS ###
###########################


# Short read alignment
def broadcast_ref_index(suffixes, reference):
    if check_existence(suffixes):
        click.echo("\n%s The following index files already exist:\n%s\nSkipping reference genome indexing..."
                   % (timestamp(), '\n'.join(suffixes)))
    else:
        click.echo("\n%s Need to generate index files for %s! Indexing reference genome %s..." % (timestamp(),
                                                                                                  reference, reference))


def broadcast_alignment(reads, reference, aligned_reads):
    if check_existence(aligned_reads):
        click.echo("\n%s Aligned reads file %s already exists!\nSkipping read alignment step..." % (timestamp(),
                                                                                                   aligned_reads))
    else:
        if len(reads) == 1:
            click.echo("\n%s Aligning read %s against the reference genome %s..." % (timestamp(), reads[0], reference))
        else:
            click.echo("\n%s Aligning the following read(s) against reference genome %s:\n%s" % (timestamp(),
                                                                                                    reference,
                                                                                                    '\n'.join(reads)))


def broadcast_sort_convert(sam_file):
    click.echo("\n%s Sorting and converting %s to the BAM format..." % (timestamp(), sam_file))


# Variant calling and post-alignment processing
def broadcast_calling(sample_list, variant_caller):
    click.echo("\n%s Calling variants with %s on the the following samples:\n%s" % (timestamp(), variant_caller,
                                                                              '\n'.join(sample_list)))


def broadcast_faidx(reference):
    faidx_file = reference + '.fai'
    if check_existence([faidx_file]):
        click.echo("\n%s Reference %s already has a faidx index file %s!\nSkipping faidx indexing." % (timestamp(),
                                                                                                       reference, faidx_file))
    else:
        click.echo("\n%s Generating faidx index %s for reference file %s..." % (timestamp(), faidx_file, reference))


def broadcast_dictionary(dict_file):
    if check_existence([dict_file]):
        click.echo("\n%s Dictionary file %s already exists!\nSkipping reference genome dictionary file generation..."
                   % (timestamp(), dict_file))
    else:
        click.echo("\n%s Generating reference genome dictionary %s..." % (timestamp(), dict_file))


def broadcast_indexing(sample_file):
    click.echo("\n%s Indexing sample file %s..." % (timestamp(), sample_file))