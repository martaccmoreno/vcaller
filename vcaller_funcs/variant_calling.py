from vcaller_funcs.auxiliary import *


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