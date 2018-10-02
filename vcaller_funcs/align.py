from vcaller_funcs.auxiliary import *


def func_align_bowtie2(output, reference: str, read1: str, read2: str = '', no_clean: bool = False) -> None:
    """
    Invokes command line tool Bowtie 2 to align read(s) against a reference sequence.
    :param output: Name of the final output file. Should end with suffix .bam
    :param reference: Path to the reference sequence against which to align the reads(s).
    :param read1: Path to a single-ended read / first read of a pair.
    :param read2: Path to the second read of the pair (Optional)
    :param no_clean: Whether or not intermediate files should be removed.
    :return: None
    """
    broadcast_step("alignment")
    suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    suffixes = [remove_suffix(reference) + suffix for suffix in suffix_list]
    broadcast_ref_index(suffixes, reference)
    if not check_existence(suffixes):
        index_args = [config['filePaths']['bowtie2'] + '/bowtie2-build', reference, remove_suffix(reference)]
        index_proc = subprocess.Popen(index_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        index_spinner = PieSpinner("\n%s Indexing reference genome %s " % (timestamp(), reference))
        progress_spinner(index_proc, index_spinner)

    sam_output = replace_suffix(output, 'sam')
    broadcast_alignment([read1, read2], reference, sam_output)
    align_args = [config['filePaths']['bowtie2'] + '/bowtie2', '-x', remove_suffix(reference), '-S', sam_output,
                  read1]
    if not check_existence([sam_output]):
        if read2 is not None:
            align_args += [read2]
        align_proc = subprocess.Popen(align_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        align_spinner = PieSpinner("Aligning ")
        progress_spinner(align_proc, align_spinner)

    if not check_existence([output]):
        sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                     os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
        broadcast_sorting(sam_output)
        try:
            subprocess.run(sort_args, check=True)
        except subprocess.CalledProcessError as e:
            print(e)

    if no_clean is False:
        cleanup([sam_output])


def func_align_bwa(output: str, reference: str, read1: str, read2: str = '', no_clean: bool = False) -> None:
    """
    Invokes command line tool BWA to align read(s) against a reference sequence.
    :param output: Name of the final output file. Should end with suffix .bam
    :param reference: Path to the reference sequence against which to align the reads(s).
    :param read1: Path to a single-ended read / first read of a pair.
    :param read2: Path to the second read of the pair (Optional).
    :param no_clean: Whether or not intermediate files should be removed.
    :return: None
    """
    suffix_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    suffixes = [remove_suffix(reference) + suffix for suffix in suffix_list]
    broadcast_ref_index(suffixes, reference)
    if not check_existence([reference + suffix for suffix in suffix_list]):
        index_args = ['bwa', 'index', reference]
        index_proc = subprocess.Popen(index_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        index_spinner = PieSpinner("%s Indexing reference genome %s " % (timestamp(), reference))
        progress_spinner(index_proc, index_spinner)

    sam_output = replace_suffix(output, 'sam')
    broadcast_alignment([read1, read2], reference, sam_output)
    align_args = ['bwa', 'mem', '-M', reference, read1]
    if check_existence([sam_output]):
        if read2 is not '':
            align_args += [read2]
        with open(sam_output, "w+") as align_out:
            align_proc = subprocess.Popen(align_args, stderr=subprocess.PIPE, stdout=align_out)
            align_spinner = PieSpinner("Aligning ")
            progress_spinner(align_proc, align_spinner)

    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                 os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
    broadcast_sorting(sam_output)
    if not check_existence([output]):
        try:
            subprocess.run(sort_args, check=True)
        except subprocess.CalledProcessError as e:
            print(e)

    if no_clean is False:
        cleanup([sam_output])


def func_align_tmap(output: str, reference: str, read1: str, read2: str = '') -> None:
    """
    Invokes command line tool TMAP to align read(s) against a reference sequence.
    :param output: Name of the final output file. Should end with suffix .bam
    :param reference: Path to the reference sequence against which to align the reads(s).
    :param read1: Path to a single-ended read / first read of a pair.
    :param read2: Path to the second read of the pair (Optional).
    :return: None
    """
    suffix_list = ['.tmap.anno', '.tmap.bwt', '.tmap.pac', '.tmap.sa']
    suffixes = [reference + suffix for suffix in suffix_list]
    broadcast_ref_index(suffixes, reference)
    if not check_existence([reference + suffix for suffix in suffix_list]):
        index_args = [config['filePaths']['tmap'], 'index', '-f', reference]
        index_proc = subprocess.Popen(index_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        index_spinner = PieSpinner("\n%s Indexing reference genome %s " % (timestamp(), reference))
        progress_spinner(index_proc, index_spinner)

    if not check_existence([output]):
        if read2 is None:
            align_args = [config['filePaths']['tmap'], 'map1', '-o', '2', '-f', reference, '-r', read1]
            if read2 is not None:
                align_args += [read2]
            if 'gz' in read1.split('.') or 'gz' in read2.split('.'):
                align_args += ['--input-gz']
            broadcast_alignment([read1, read2], reference, output)
            if not check_existence([output]):
                with open(output, "w+") as align_out:
                    align_proc = subprocess.Popen(align_args, stderr=subprocess.PIPE, stdout=align_out)
                    align_spinner = PieSpinner("Aligning ")
                    progress_spinner(align_proc, align_spinner)