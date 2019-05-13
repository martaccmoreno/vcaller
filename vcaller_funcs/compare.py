from vcaller_funcs.auxiliary import *

def func_compare(output_dir, bed_file, evaluation_regions, score_field, sample, reference, baseline, calls, no_clean):
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
        tabix_index([baseline])
        baseline += '.gz'
    if not check_existence(calls + '.tbi'):
        tabix_index([calls])
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