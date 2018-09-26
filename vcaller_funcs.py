import json
from auxiliary_funcs import *


#####################
### IMPORT CONFIG ###
#####################

current_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(current_dir, 'config.json'), 'r') as data_file:
    config = json.load(data_file)


######################
### MAIN FUNCTIONS ###
######################

### ALIGNERS ####
def func_align_bowtie2(output, reference, read1, read2='', no_clean=False):
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

    sort_args = ['samtools', 'sort', '-O', 'bam', '-o', output, '-T',
                 os.path.join('/tmp/', replace_suffix(os.path.basename(output), 'tmp')), sam_output]
    broadcast_sort_convert(sam_output)
    if not check_existence([output]):
        try:
            subprocess.run(sort_args)
        except subprocess.CalledProcessError as e:
            print(e)

    if no_clean is False:
        cleanup(sam_output)


def func_align_bwa(output, reference, read1, read2='', no_clean=False):
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
    if not check_existence([sam_output]):
        subprocess.run(sort_args)

    if no_clean is False:
        cleanup(sam_output)


def func_align_tmap(output, reference, read1, read2=''):
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


### VARIANT CALLERS ###
def func_call_bcftools(output, exome_regions, reference, sample1, sample2='', count_orphans=False, no_clean=False):
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


def func_call_freebayes(output, exome_regions, reference, sample1, sample2=''):
    sample_list = [sample1] + [s for s in sample2]

    broadcast_calling(sample_list, "Freebayes")
    call_args = [config['filePaths']['freebayes'], '-f', reference] + sample_list
    with open(output, 'w+') as call_out:
        subprocess.run(call_args, stdout=call_out)

    if exome_regions:
        bed_intersect(output, exome_regions, clean=True)


def func_call_gatk(output, dbsnp, exome_regions, reference, sample1, sample2=''):
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
        if check_existence([replace_suffix(smpl, 'bai')]) or check_existence(smpl + '.bai'):
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


def func_call_tvc(output_dir, exome_regions, reference, sample):
    broadcast_calling([sample], "TVC")
    call_args = [config['filePaths']['tvc'], '-i', sample, '-r', reference, '-o', output_dir]
    if exome_regions:
        call_args += ['-b', exome_regions]
    subprocess.run(call_args)


def func_call_varscan2(output, count_orphans, exome_regions, no_clean, reference, sample1, sample2=''):
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


### POST-ALIGNMENT PROCESSING ####
def func_process(output_name, output_dir, readgroup_info, add_known_indels, known_indels, known_snps, reference,
                 sample, no_clean=False):
    if output_dir is None:
        output_dir = os.path.dirname(output_name)
    smpl_name, smpl_extension = remove_suffix(os.path.basename(sample)), '.' + sample.split('.')[-1]

    if not check_existence([reference + '.fai']):
        broadcast_faidx(reference)
        subprocess.run(['samtools', 'faidx', reference])
    dict_file = remove_suffix(reference) + '.dict'
    if not check_existence([dict_file]):
        broadcast_dictionary(dict_file)
        subprocess.run([config['filePaths']['gatk4'], 'CreateSequenceDictionary', '-R', reference])

    if (smpl_extension is not '.bam') or (subprocess.getstatusoutput('samtools index ' + sample)[0] != 0):
        sorted_output = replace_suffix(sample, 'bam')
        if not check_existence([sorted_output]):
            broadcast_sort_convert(sample)
            sort_args = ['samtools', 'sort', '-O', 'bam', '-o', sorted_output, '-T',
                         os.path.join('/tmp/', smpl_name + '.temp'), sample]
            subprocess.run(sort_args)
            broadcast_indexing(sorted_output)
            subprocess.run(['samtools', 'index', sorted_output])
        sample = sorted_output

    if check_existence(sample + '.bai'):  # to avoid unnecessary 'file not found' errors
        subprocess.run(['rm', sample + '.bai'])

    # dedupping
    dup_output = os.path.join(output_dir, smpl_name + '.DUP.bam')
    dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', sample,
                '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                '-M', os.path.join(output_dir, smpl_name + '.metrics')]
    if readgroup_info:
        # More info on RGs: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        rg_output = os.path.join(output_dir, smpl_name + '.RG.bam')
        if not check_existence(rg_output):
            click.echo("\n%s Adding Read Group information to %s..." % (timestamp(), sample))
            rg_info = readgroup_info.split(',')
            read_groups = {'ID': rg_info[0].split(':')[1], 'PU': rg_info[1].split(':')[1],
                           'PL': rg_info[2].split(':')[1], 'SM': rg_info[3].split(':')[1],
                           'LB': rg_info[4].split(':')[1]}
            rg_args = [config['filePaths']['gatk4'], 'AddOrReplaceReadGroups', '-I', sample,
                       '-O', rg_output, '-RGID', read_groups['ID'], '-RGLB', read_groups['LB'],
                       '-RGPL', read_groups['PL'].upper(), '-RGPU', read_groups['PU'], '-RGSM', read_groups['SM'],
                       '-SO', 'coordinate']
            subprocess.run(rg_args)
        dup_output = replace_suffix(rg_output, '.DUP') + '.bam'
        dup_args = [config['filePaths']['gatk4'], 'MarkDuplicates', '-I', rg_output,
                    '-O', dup_output, '-REMOVE_DUPLICATES', 'True',
                    '-M', os.path.join(output_dir, smpl_name + '.metrics')]

    if not check_existence(dup_output):
        click.echo("\n%s Marking and removing duplicates for %s..." % (timestamp(), sample))
        subprocess.run(dup_args)

    if not check_existence(dup_output+'.bai'):
        broadcast_indexing(dup_output)
        subprocess.run(['samtools', 'index', dup_output])

    # Realign around indels; using gatk3 because of this step
    # https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
    intervals_output = os.path.join(output_dir, smpl_name + '.intervals')
    if not check_existence([intervals_output]):
        click.echo("\n%s Creating indel realignment intervals for %s..." % (timestamp(), dup_output))
        intervals_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'RealignerTargetCreator', '-R', reference,
                          '-I', dup_output, '-o', intervals_output, '--known',
                          known_indels]
        subprocess.run(intervals_args)
    realign_output = replace_suffix(dup_output, 'RLGN') + '.bam'
    if not check_existence([realign_output]):
        click.echo("\n%s Applying indel realignment based on the intervals for %s..." % (timestamp(), dup_output))
        realign_args = ['java', '-jar', config['filePaths']['gatk3'], '-T', 'IndelRealigner', '-R', reference,
                        '-I', dup_output, '-targetIntervals', intervals_output, '-known', known_indels] + \
                       flatten_list([['-known'] + [add_known_indels[i]] for i in range(len(add_known_indels))]) + \
                       ['-o', realign_output]
        subprocess.run(realign_args)

    # BQSR
    table_output = os.path.join(output_dir, smpl_name + '.table')
    if not check_existence([table_output]):
        click.echo("\n%s Creating base score recalibration table for %s..." % (timestamp(), realign_output))
        table_args = [config['filePaths']['gatk4'], 'BaseRecalibrator', '-R', reference,
                      '--known-sites', known_snps] + ['-I', realign_output, '-O', table_output]
        subprocess.run(table_args)
    if output_name is None:
        bqsr_output = os.path.join(output_dir, smpl_name + '.processed.bam')
    else:
        bqsr_output = output_name
    if not check_existence([bqsr_output]):
        click.echo("\n%s Running base score recalibration on %s..." % (timestamp(), realign_output))
        bqsr_args = [config['filePaths']['gatk4'], 'ApplyBQSR',
                     '-I', realign_output, '-bqsr', table_output, '-O', bqsr_output]
        subprocess.run(bqsr_args)

    # clean up intermediary files -- but only after we have the final file
    # still not working 100% -- IT STILL RUNS EVEN IF BSQR WASN'T SUCCESSFUL
    if check_existence([bqsr_output]) and no_clean is False:
        files_to_rmv = [dup_output, intervals_output, realign_output, table_output,
                        os.path.join(output_dir, smpl_name + '.metrics')]
        if readgroup_info:
            files_to_rmv += [rg_output]
        index_files = [replace_suffix(item, 'bai') for item in [file for file in files_to_rmv if '.bam' in file]]
        cleanup(files_to_rmv+index_files)


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