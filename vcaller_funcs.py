from auxiliary_funcs import *



#####################
### IMPORT CONFIG ###
#####################

config = import_config()

######################
### MAIN FUNCTIONS ###
######################

### ALIGNERS ####
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
        subprocess.run(index_args)

    sam_output = replace_suffix(output, 'sam')
    broadcast_alignment([read1, read2], reference, sam_output)
    align_args = [config['filePaths']['bowtie2'] + '/bowtie2', '-x', remove_suffix(reference), '-S', sam_output,
                  read1]
    if not check_existence([sam_output]):
        if read2 is not None:
            align_args += [read2]
        subprocess.run(align_args)

    if not check_existence([output]):
        sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                     os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
        broadcast_sort_convert(sam_output)
        try:
            subprocess.run(sort_args)
        except subprocess.CalledProcessError as e:
            print(e)

    if no_clean is False:
        cleanup(sam_output)


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
        subprocess.run(['bwa', 'index', reference])

    sam_output = replace_suffix(output, 'sam')
    broadcast_alignment([read1, read2], reference, sam_output)
    align_args = ['bwa', 'mem', '-M', reference, read1]
    if check_existence([sam_output]):
        if read2 is not '':
            align_args += [read2]
        with open(sam_output, "w+") as align_out:
            subprocess.run(align_args, stdout=align_out)

    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                 os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
    broadcast_sort_convert(sam_output)
    if not check_existence([output]):
        try:
            subprocess.run(sort_args, check=True)
        except subprocess.CalledProcessError as e:
            print(e)

    if no_clean is False:
        cleanup(sam_output)


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
        subprocess.run(index_args)

    if not check_existence([output]):
        if read2 is None:  # if read is single-ended
            align_args = [config['filePaths']['tmap'], 'map1', '-o', '2', '-f', reference, '-r', read1]
            if read2 is not None:
                align_args += [read2]
            if 'gz' in read1.split('.') or 'gz' in read2.split('.'):
                align_args += ['--input-gz']
            broadcast_alignment([read1, read2], reference, output)
            if not check_existence([output]):
                with open(output, "w+") as align_out:
                    subprocess.run(align_args, stdout=align_out)


### POST-ALIGNMENT PROCESSING ####
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

### VARIANT CALLERS ###
def func_call_bcftools(output: str, exome_regions: str, reference: str, sample1: str, sample2: str = '',
                       count_orphans: bool = False, no_clean: bool = False) -> None:
    """
    Invokes BCFtools to call variants on aligned reads (samples).
    :param output: Name of the final output file. Should end with suffix .vcf
    :param exome_regions: Bed file defining exome regions from which to extract calls for the final output file.
    :param reference: Reference sequence against which to compare aligned reads.
    :param sample1: A data set containing aligned reads.
    :param sample2: Another data set contained aligned reads (Optional; needed for multi-sample/join variant calling)
    :param count_orphans: Whether or not to skip anomalous read pairs in variant calling.
    :param no_clean: Whether or not intermediate files should be removed.
    :return: None
    """
    sample_list = [sample1] + [s for s in sample2]
    bcf_output = replace_suffix(output, 'bcf')

    click.echo("Calculating genotype likelihoods for %s..." % ', '.join(sample_list))
    if count_orphans:
        mpileup_args = ['bcftools', 'mpileup', '-AOb', '-o', bcf_output, '-f', reference] + sample_list
    else:
        mpileup_args = ['bcftools', 'mpileup', '-Ob', '-o', bcf_output, '-f', reference] + sample_list
    subprocess.run(mpileup_args)

    broadcast_calling(sample_list, "BCFtools")
    call_args = ['bcftools', 'call', '-vmO', 'v', '-o', output, bcf_output]
    subprocess.run(call_args)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)
    if no_clean is False:
        cleanup(bcf_output)


def func_call_freebayes(output: str, exome_regions: str, reference: str, sample1: str, sample2: str = '') -> None:
    """
    Invokes Freebayes to call variants on aligned reads (samples).
    :param output: Name of the final output file. Should end with suffix .vcf
    :param exome_regions: Bed file defining exome regions from which to extract calls for the final output file.
    :param reference: Reference sequence against which to compare aligned reads.
    :param sample1: A data set containing aligned reads.
    :param sample2: Another data set contained aligned reads (Optional; needed for multi-sample/join variant calling)
    :return: None
    """
    sample_list = [sample1] + [s for s in sample2]

    broadcast_calling(sample_list, "Freebayes")
    call_args = [config['filePaths']['freebayes'], '-f', reference] + sample_list
    with open(output, 'w+') as call_out:
        subprocess.run(call_args, stdout=call_out)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)


def func_call_haplotypecaller(output: str, dbsnp: str, exome_regions: str, reference: str, sample1: str,
                              sample2: str = '') -> None:
    """
    Invokes GATK's HaplotypeCaller to call variants on aligned reads (samples).
    :param output: Name of the final output file. Should end with suffix .vcf
    :param dbsnp: Path to a dbsnp file containing a list of ids for known SNPs.
    :param exome_regions: Bed file defining exome regions from which to extract calls for the final output file.
    :param reference: Reference sequence against which to compare aligned reads.
    :param sample1: A data set containing aligned reads.
    :param sample2: Another data set contained aligned reads (Optional; needed for multi-sample/join variant calling)
    :return: None
    """
    broadcast_faidx(reference)
    if not check_existence([reference + '.fai']):
        subprocess.run(['samtools', 'faidx', reference])

    dict_file = replace_suffix(reference, 'dict')
    dict_vars = ['java', '-jar', config['filePaths']['picard'], 'CreateSequenceDictionary', 'R=%s' % reference,
                 'O=%s' % dict_file]
    broadcast_dictionary(dict_file)
    if not check_existence([dict_file]):
        subprocess.run(dict_vars)

    sample_list = [sample1] + [s for s in sample2]
    for smpl in sample_list:
        if check_existence([replace_suffix(smpl, 'bai')]) or check_existence([smpl + '.bai']):
            click.echo("\n%s The sample index .bai file already exists!\nSkipping sample indexing..." % timestamp())
        else:
            click.echo('\n%s Need to generate sample index .bai file!\nIndexing sample file %s...' % (timestamp(), smpl))
            subprocess.run(['samtools', 'index', smpl])

    broadcast_calling(sample_list, "HaplotypeCaller")
    if dbsnp is None:
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference] + \
                    flatten_list([['-I'] + [sample_list[i]] for i in range(len(sample_list))]) + ['-O', output]
    else:
        gatk_args = [config['filePaths']['gatk4'], 'HaplotypeCaller', '-R', reference] + \
                    flatten_list([['-I'] + [sample_list[i]] for i in range(len(sample_list))]) + \
                    ['--dbsnp', dbsnp, '-O', output]
    subprocess.run(gatk_args)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)


def func_call_tvc(output_dir: str, exome_regions: str, reference: str, sample: str) -> None:
    """
    Invokes TVC to call variants on aligned reads (sample).
    :param output_dir: Name of the output directory.
    :param exome_regions: Bed file defining exome regions from which to extract calls for the final output file.
    :param reference: Reference sequence against which to compare aligned reads.
    :param sample: A data set containing aligned reads.
    :return: None
    """
    broadcast_calling([sample], "TVC")
    call_args = [config['filePaths']['tvc'], '-i', sample, '-r', reference, '-o', output_dir]
    if exome_regions:
        call_args += ['-b', exome_regions]
    subprocess.run(call_args)


def func_call_varscan2(output: str, exome_regions: str, reference: str, sample1: str,
                       sample2: str = '', count_orphans: bool = False, no_clean: bool = True) -> None:
    """
    Invokes VarScan 2 to call variants on aligned reads (samples).
    :param output: Name of the final output file. Should end with suffix .vcf
    :param exome_regions: Bed file defining exome regions from which to extract calls for the final output file.
    :param reference: Reference sequence against which to compare aligned reads.
    :param sample1: A data set containing aligned reads.
    :param sample2: Another data set contained aligned reads (Optional; needed for multi-sample/join variant calling)
    :param count_orphans: Whether or not to skip anomalous read pairs in variant calling.
    :param no_clean: Whether or not intermediate files should be removed.
    :return: None
    """
    sample_list = [sample1] + [s for s in sample2]
    click.echo("\n%s Creating mpipleup file for the following samples:\n%s..." % (timestamp(), '\n'.join(sample_list)))
    mpileup_file = replace_suffix(output, 'pileup')
    if count_orphans:
        pileup_args = ['samtools', 'mpileup', '-A', '-f', reference] + sample_list
    else:
        pileup_args = ['samtools', 'mpileup', '-f', reference] + sample_list
    with open(mpileup_file, 'w+') as pileup_out:
        subprocess.run(pileup_args, stdout=pileup_out)

    broadcast_calling(sample_list, "VarScan 2")
    call_args = ['java', '-jar', config['filePaths']['varscan2'], 'mpileup2cns', mpileup_file, '--output-vcf', '1',
                 '--variants', '1', '--p-value', '0.10', '--min-coverage', '2']
    with open(output, 'w+') as call_out:
        subprocess.run(call_args, stdout=call_out)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)
    if no_clean is False:
        cleanup(mpileup_file)


### VARIANT CALLING BENCHMARKING ###
def func_compare(output_dir, bed_file, evaluation_regions, score_field, sample, no_clean, reference, baseline, calls):
    # The reference genome must be converted to SDF
    sdf_ref = os.path.join(replace_suffix(reference, 'sdf'))
    if os.path.isdir(sdf_ref) is False:
        click.echo("\n%s Converting reference genome %s to the SDF format..." % (timestamp(), reference))
        fastq2sdf_args = [config['filePaths']['rtg'], 'format', '-o', sdf_ref, reference]
        subprocess.run(fastq2sdf_args)
    else:
        click.echo("\n%s Reference genome %s has already been convert to the SDF format as %s" %
                   (timestamp(), os.path.basename(reference), sdf_ref))

    # Check if baseline and calls are tabix-indexed as this is a must for using vcfeval
    if not check_existence(baseline + '.tbi'):
        tabix_index(baseline)
        baseline += '.gz'
    if not check_existence(calls + '.tbi'):
        tabix_index(calls)
        calls += '.gz'

    click.echo("\n%s Creating directory %s..." % (timestamp(), output_dir))

    # Create GA4GH-compliant annotated VCFs
    rtg_out = os.path.join(output_dir, os.path.split(output_dir)[-1] + '-vcfeval')
    if not os.path.isdir(rtg_out):
        click.echo("\n%s Comparing baseline %s against call set %s using vcfeval..." % (timestamp(), baseline, calls))
        rtg_args = [config['filePaths']['rtg'], 'vcfeval', '-o', rtg_out, '--vcf-score-field',
                    score_field, '--template', sdf_ref, '--baseline', baseline, '--calls', calls, '-m', 'ga4gh']
        if bed_file is not None: rtg_args += ['--bed-regions', bed_file]
        if evaluation_regions is not None: rtg_args += ['--evaluation-regions', evaluation_regions]
        if sample is not None: rtg_args += ['--sample', sample]
        subprocess.run(rtg_args)

    click.echo("\n%s Moving to directory %s..." % (timestamp(), output_dir))
    initial_path = os.getcwd()
    os.chdir(output_dir)
    click.echo("\n%s Running qfy.py on %s..." % (timestamp(), os.path.join(rtg_out, 'output.vcf.gz')))
    qfy_args = [config['filePaths']['qfy.py'], '-t', 'ga4gh', '--verbose', '--adjust-conf-regions',
                os.path.normpath(os.path.join(initial_path, evaluation_regions)), '--reference',
                os.path.normpath(os.path.join(initial_path, reference)), '-o', os.path.split(output_dir)[-1],
                '--write-vcf', '--write-counts',
                os.path.join(os.path.basename(rtg_out), 'output.vcf.gz')]
    subprocess.run(qfy_args)
    click.echo("\n%s Returning to %s..." % (timestamp(), initial_path))
    os.chdir(initial_path)

    if no_clean is False:
        cleanup(rtg_out)