import json
import click
import os
import subprocess
import time
from progress.spinner import PieSpinner # new dependency
from progress.bar import IncrementalBar
import re


###########################
### AUXILIARY FUNCTIONS ###
###########################
def import_config(config_json: str) -> dict:
    """
    Imports a json configuration file as a dictionary.
    :param config_json: The json file containing configuration info.
    :return: A dictionary encoding the information provided by the json file.
    """
    basepath = os.path.dirname(__file__)
    filepath = os.path.abspath(os.path.join(basepath, "..",  config_json))
    with open(filepath, 'r') as json_file:
        return json.load(json_file)


def flatten_list(list_of_list: list) -> list:
    """Flatten a list of list into a single list."""
    return [item for sublist in list_of_list for item in sublist]


def remove_suffix(file_name: str) -> str:
    """Remove the suffix of a filename, e.g. 'reference.fa' becomes 'reference'"""
    return '.'.join(file_name.split('.')[:-1])


def replace_suffix(filename: str, new_suffix: str) -> str:
    """
    Replaces a filename's suffix with another user-specified suffix.
    '"""
    if new_suffix[0] == '.':
        return remove_suffix(filename) + new_suffix
    else:
        return remove_suffix(filename) + '.' + new_suffix


def check_existence(filename_list: list) -> bool:
    """
    Check if files with the filenames in the list already exist.
    """
    if sum([os.path.isfile(ifile) for ifile in filename_list]) == len(filename_list):
        return True
    else:
        return False


def tabix_index(bgzipped_files: list) -> None:
    """
    Tabix index bgzipped files, or bgzip then index regular files, or gunzip then bgzip then index gzipped files.
    :param bgzipped_files: list of files to index.
    :return: None
    """
    for file in bgzipped_files:
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


def cleanup(files_or_dirs: list) -> None:
    """
    Permanently remove one or more files or directories.
    :param files_or_dirs: a list of files and/or directories to remove
    :return: None
    """
    click.echo('\nRemoving the following files/directories: %s' % ', '.join(files_or_dirs))
    for file_or_dir in files_or_dirs:
        if os.path.isfile(file_or_dir):
            subprocess.run(['rm', file_or_dir])
        elif os.path.isdir(file_or_dir):
            subprocess.run(['rm', '-r', file_or_dir])
        else:
            click.echo('\nInvalid input: %s is neither a file nor a directory...' % file_or_dir)


def bed_intersect(vcf: str, bed: str, out: str = None, clean: bool = False) -> None:
    """
    Intersect a vcf file with a bed file, obtaining a second vcf file with only the regions defined in the bed file.
    """
    if out is None:
        out = remove_suffix(vcf)+'_exome.vcf'
    with open(out, 'w+') as out:
        click.echo('\nIntersecting vcf %s with bed file regions in %s...' % (vcf, bed))
        subprocess.run(['bedtools', 'intersect', '-header', '-a', vcf, '-b', bed], stdout=out)
    if clean:
        cleanup([vcf])


def timestamp() -> str:
    """
    Prints a pretty timestamp.
    """
    return time.strftime("[%Y-%m-%d %H:%M]")


###########################
### BROADCAST FUNCTIONS ###
###########################
def broadcast_step(step):
    click.echo("%s *** Starting %s step. ***" % (timestamp(), step))


def broadcast_error(error_code, command, error_message):
    click.echo("%s An error has been found during Vcaller's execution:\n" % timestamp())
    if error_message:
        click.echo(error_message)
    raise subprocess.CalledProcessError(error_code, command)


def broadcast_ref_index(suffixes: list, reference: str):
    if check_existence(suffixes):
        click.echo("%s The following index files already exist:\n%s" % (timestamp(), '\n'.join(suffixes)))
        click.echo("%s Skipping reference genome indexing." % timestamp())
    else:
        click.echo("%s Need to generate index files for %s!" % (timestamp(), reference))


def broadcast_alignment(reads, reference, aligned_reads):
    if check_existence(aligned_reads):
        click.echo("\n%s Aligned reads file %s already exists!" % (timestamp(), aligned_reads))
        click.echo("\n%s Skipping read alignment step." % timestamp())
    else:
        if len(list(filter(None, reads))) == 1:
            click.echo("%s Aligning read %s against the reference genome %s..." % (timestamp(), reads[0], reference))
        else:
            click.echo("%s Aligning the following read(s) against reference genome %s:\n%s" % (timestamp(), reference,
                                                                                                 '\n'.join(reads)))


def broadcast_sorting(unsorted_file):
    click.echo("%s Sorting and converting %s to the BAM format." % (timestamp(), unsorted_file))


def broadcast_calling(sample_list, variant_caller):
    click.echo("%s Calling variants with %s on the the following samples:\n%s" % (timestamp(), variant_caller,
                                                                              '\n'.join(sample_list)))


def broadcast_faidx(reference):
    faidx_file = reference + '.fai'
    if check_existence([faidx_file]):
        click.echo("%s Reference %s already has a faidx index file %s!\nSkipping faidx indexing." % (timestamp(),
                                                                                                       reference, faidx_file))
    else:
        click.echo("%s Generating faidx index %s for reference file %s..." % (timestamp(), faidx_file, reference))


def broadcast_dictionary(dict_file):
    if check_existence([dict_file]):
        click.echo("%s Dictionary file %s already exists!\nSkipping reference genome dictionary file generation..."
                   % (timestamp(), dict_file))
    else:
        click.echo("%s Generating reference genome dictionary %s..." % (timestamp(), dict_file))


def broadcast_indexing(sample_file):
    click.echo("%s Indexing sample file %s." % (timestamp(), sample_file))


############################
### PROGRESS MEASUREMENT ###
############################
def progress_spinner(process: subprocess.Popen, spinner: PieSpinner):
    """
    Create a spinner that is updated as a given process runs. Raise the exit status if any error is found during
    the process's execution.
    :param process: the process's Popen object
    :param spinner: A Spinner object to update
    :return: A progress message that updates as the process runs.
    """
    error_message = ""
    spinner.start()
    spinner.next()
    while True:
        line = process.stderr.readline().decode('utf-8')
        error_message += line
        if not line:
            break
        else:
            spinner.next()
    spinner.finish()
    process.communicate()
    if process.returncode > 0:
        broadcast_error(process.returncode, process, error_message)


def progress_bar(process: subprocess.Popen, bar: IncrementalBar):
    """
    Create a loading bar that is updated as a given process runs based on its completion %.
    Raise the exit status if any error is found during the process's execution.
    :param process: the process's Popen object
    :param bar: A Bar object to update
    :return: A progress message that updates as the process runs.
    """
    error_message = ""
    bar.start()
    bar.next()
    while True:
        line = process.stderr.readline().decode('utf-8')
        error_message += line
        if not line:
            break
        try:
            percent_match = re.search("[1-9]+\.[1-9]%", line)
            if percent_match:
                percent = float(percent_match.group(0).replace("%", ""))
                bar.goto(percent)
        except IndexError:
            pass
    bar.goto(100)
    bar.finish()
    process.communicate()
    if process.returncode > 0:
        broadcast_error(process.returncode, process, error_message)


#########################
### IMPORT TOOL PATHS ###
#########################
config = import_config('tool_paths.json')