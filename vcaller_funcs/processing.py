from vcaller_funcs.auxiliary import *
import subprocess


def func_process(output_name: str, output_dir: str, readgroup_info: str, add_known_indels: str, known_indels: str,
                 known_snps: str, reference: str, sample: str, no_clean: bool = False) -> None:
    """
    Perform post-alignment steps on a data set containing aligned reads.
    :param output_name:
    :param output_dir:
    :param readgroup_info:
    :param add_known_indels:
    :param known_indels:
    :param known_snps:
    :param reference:
    :param sample:
    :param no_clean:
    :return:
    """
    broadcast_step("post-alignment processing")
    if output_dir is None:
        output_dir = os.path.dirname(output_name)
    smpl_name, smpl_extension = remove_suffix(os.path.basename(sample)), '.' + sample.split('.')[-1]

    if not check_existence([reference + '.fai']):
        try:
            broadcast_faidx(reference)
            subprocess.run(['samtools', 'faidx', reference], check=True) # no need to capture std
        except subprocess.CalledProcessError as e:
            print(e)
    dict_file = remove_suffix(reference) + '.dict'
    if not check_existence([dict_file]):
        try:
            broadcast_dictionary(dict_file)
            subprocess.run([config['filePaths']['gatk4'], 'CreateSequenceDictionary', '-R', reference], check=True)
        except subprocess.CalledProcessError as e:
            click.echo(e)

    if (smpl_extension is not '.bam') or (subprocess.getstatusoutput('samtools index ' + sample)[0] != 0):
        sorted_output = replace_suffix(sample, 'bam')
        if not check_existence([sorted_output]):
            broadcast_sort_convert(sample)
            try:
                sort_args = ['samtools', 'sort', '-O', 'bam', '-o', sorted_output, '-T',
                             os.path.join('/tmp/', smpl_name + '.temp'), sample]
                subprocess.run(sort_args, check=True)
            except subprocess.CalledProcessError as e:
                click.echo(e)
            try:
                broadcast_indexing(sorted_output)
                subprocess.run(['samtools', 'index', sorted_output])
            except subprocess.CalledProcessError as e:
                click.echo(e)
        sample = sorted_output

    # if check_existence(sample + '.bai'):  # to avoid unnecessary 'file not found' errors
    #     subprocess.run(['rm', sample + '.bai'])

    # dedupping
    dup_output = os.path.join(output_dir, smpl_name + '.DUP.bam')
    dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', sample,
                '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                '-M', os.path.join(output_dir, smpl_name + '.metrics')]
    if not check_existence(dup_output):
        dedup_proc = subprocess.Popen(dup_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        dedup_spinner = Spinner("\n%s Marking and removing duplicates for %s " % (timestamp(), sample))
        progress_spinner(dedup_proc, dedup_spinner)

    # More info on RGs: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    rg_output = replace_suffix(dup_output, 'RG') + '.bam'
    if not check_existence(rg_output):
        rg_info = readgroup_info.split(',')
        read_groups = {'ID': rg_info[0].split(':')[1], 'PU': rg_info[1].split(':')[1],
                       'PL': rg_info[2].split(':')[1], 'SM': rg_info[3].split(':')[1],
                       'LB': rg_info[4].split(':')[1]}
        rg_args = [config['filePaths']['gatk4'], 'AddOrReplaceReadGroups', '-I', dup_output,
                   '-O', rg_output, '-RGID', read_groups['ID'], '-RGLB', read_groups['LB'],
                   '-RGPL', read_groups['PL'].upper(), '-RGPU', read_groups['PU'], '-RGSM', read_groups['SM'],
                   '-SO', 'coordinate']
        rg_proc = subprocess.Popen(rg_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        rg_spinner = Spinner("\n%s Adding Read Group information to %s " % (timestamp(), sample))
        progress_spinner(rg_proc, rg_spinner)
        if not check_existence(dup_output+'.bai'):
            broadcast_indexing(rg_output)
            subprocess.run(['samtools', 'index', rg_output])
    # elif not check_existence(dup_output+'.bai'):
    #     broadcast_indexing(dup_output)
    #     subprocess.run(['samtools', 'index', dup_output])


    # Realign around indels; using gatk3 because of this step
    # https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
    # Generate intervals in which to realign around indels
    intervals_output = os.path.join(output_dir, smpl_name + '.intervals')
    intervals_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'RealignerTargetCreator', '-R', reference,
                      '-I', rg_output, '-o', intervals_output, '--known',
                      known_indels]
    if not check_existence([intervals_output]):
        intervals_proc = subprocess.Popen(intervals_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        click.echo("") # create only ONE newline, not one every time the bar updates
        intervals_bar = IncrementalBar("%s Creating indel realignment intervals for %s " % (timestamp(), rg_output),
                                       suffix='%(percent).1f%% - %(elapsed)ds')
        intervals_bar.start()
        progress_bar(intervals_proc, intervals_bar)
    # Applying indel realignment in the defined regions
    realign_output = replace_suffix(rg_output, 'RLGN') + '.bam'
    if not check_existence([realign_output]):
        realign_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'IndelRealigner', '-R', reference,
                        '-I', rg_output, '-targetIntervals', intervals_output, '-known', known_indels] + \
                       flatten_list([['-known'] + [add_known_indels[i]] for i in range(len(add_known_indels))]) + \
                       ['-o', realign_output]
        click.echo("") # create only ONE newline, not one every time the bar updates
        realign_proc = subprocess.Popen(realign_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        realign_bar = IncrementalBar("%s Applying indel realignment to %s " % (timestamp(), rg_output),
                                     suffix='%(percent).1f%% - %(elapsed)ds')
        realign_bar.start()
        progress_bar(realign_proc, realign_bar)

    # BQSR ### ADD SPINNER TO TABLE PART AND SEE WHAT THE OTHER IS LIKE ###
    table_output = os.path.join(output_dir, smpl_name + '.table')
    if not check_existence([table_output]):
        table_args = [config['filePaths']['gatk4'], 'BaseRecalibrator', '-R', reference,
                      '--known-sites', known_snps] + ['-I', realign_output, '-O', table_output]
        table_proc = subprocess.Popen(table_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        table_spinner = Spinner("\n%s Creating base score recalibration table for %s " % (timestamp(), sample))
        progress_spinner(table_proc, table_spinner)
    if output_name is None:
        bqsr_output = os.path.join(output_dir, smpl_name + '.processed.bam')
    else:
        bqsr_output = output_name
    if not check_existence([bqsr_output]):
        bqsr_args = [config['filePaths']['gatk4'], 'ApplyBQSR',
                     '-I', realign_output, '-bqsr', table_output, '-O', bqsr_output]
        bqsr_proc = subprocess.Popen(bqsr_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        bqsr_spinner = Spinner("\n%s Running base score recalibration on %s " % (timestamp(), realign_output))
        progress_spinner(bqsr_proc, bqsr_spinner)

    if check_existence([bqsr_output]) and no_clean is False:
        files_to_rmv = [dup_output, intervals_output, realign_output, table_output,
                        os.path.join(output_dir, smpl_name + '.metrics')]
        if readgroup_info:
            files_to_rmv += [rg_output]
        index_files = [replace_suffix(item, 'bai') for item in [file for file in files_to_rmv if '.bam' in file]]
        cleanup(files_to_rmv+index_files)