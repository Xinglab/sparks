import snakemake.utils

snakemake.utils.min_version('6.5.3')

configfile: 'snakemake_config.yaml'

onsuccess:
    print('workflow success')

onerror:
    print('workflow error')

DEFAULT_MEM_MB=4 * 1024  # 4 GB
DEFAULT_TIME_HOURS=12

USED_EVENT_TYPES = ['SE', 'A3SS', 'A5SS', 'RI']
ALL_EVENT_TYPES = ['SE', 'A3SS', 'A5SS', 'RI', 'MXE']

localrules: all
rule all:
    input:
        rds='{}.SPARKS.rds'.format(config['cancer_type']),


def reference_file_wildcard_constraints():
    reference_files = config.get('reference_files')
    if reference_files:
        file_names = '|'.join([re.escape(file_name)
                               for file_name in reference_files])
        without_gz = '|'.join([re.escape(file_name[:-3])
                               for file_name in reference_files
                               if file_name.endswith('.gz')])
    else:
        no_match = '^$'  # only matches empty string
        file_names = no_match
        without_gz = no_match

    return {'file_names': file_names, 'without_gz': without_gz}


def get_url_for_download_reference_file(wildcards):
    file_name = wildcards.file_name
    return config['reference_files'][file_name]['url']


rule download_reference_file:
    output:
        ref_file=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'download_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'download_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['file_names']
    params:
        url=get_url_for_download_reference_file,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl -L \'{params.url}\''
        ' -o {output.ref_file}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule unzip_reference_file:
    input:
        gz=os.path.join('references', '{file_name}.gz'),
    output:
        un_gz=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['without_gz']
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        ' gunzip -c {input.gz}'
        ' 1> {output.un_gz}'
        ' 2> {log.err}'

rule extract_transcriptome_idx:
    input:
        gz=os.path.join('references', 'homo_sapiens.tar.gz'),
    output:
        idx=os.path.join('references', 'homo_sapiens', 'transcriptome.idx'),
        to_genes=os.path.join('references', 'homo_sapiens',
                              'transcripts_to_genes.txt'),
    log:
        out=os.path.join('references', 'extract_transcriptome_idx_log.out'),
        err=os.path.join('references', 'extract_transcriptome_idx_log.err'),
    params:
        # cd to the reference dir and use relative paths
        ref_dir='references',
        gz='homo_sapiens.tar.gz',
        idx=os.path.join('homo_sapiens', 'transcriptome.idx'),
        to_genes=os.path.join('homo_sapiens', 'transcripts_to_genes.txt'),
        log_out='extract_transcriptome_idx_log.out',
        log_err='extract_transcriptome_idx_log.err',
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        ' cd {params.ref_dir}'
        ' && tar -xvf {params.gz}'
        ' {params.idx}'
        ' {params.to_genes}'
        ' 1> {params.log_out}'
        ' 2> {params.log_err}'

rule create_STAR_reference:
    input:
        gtf=os.path.join('references', config['gtf_name']),
        fasta=os.path.join('references', config['fasta_name']),
    output:
        ref=directory(os.path.join('references', 'STAR')),
    log:
        out=os.path.join('references', 'create_STAR_reference_log.out'),
        err=os.path.join('references', 'create_STAR_reference_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
    threads: config['create_STAR_reference_threads']
    resources:
        mem_mb=config['create_STAR_reference_mem_gb'] * 1024,
        time_hours=config['create_STAR_reference_time_hr'],
    shell:
        'mkdir -p {output.ref}'
        ' && {params.conda_wrapper} STAR'
        ' --runThreadN {threads}'
        ' --runMode genomeGenerate'
        ' --genomeDir {output.ref}'
        ' --genomeFastaFiles {input.fasta}'
        ' --sjdbGTFfile {input.gtf}'
        ' --sjdbOverhang 99'
        ' 1> {log.out}'
        ' 2> {log.err}'


def get_fastqs_for_sample_name(sample_name):
    for group_key in ['sample_group_1', 'sample_group_2']:
        for name, fastqs in config[group_key].items():
            if name == sample_name:
                return fastqs

    raise Exception('could not find fastqs for sample: {}'.format(sample_name))


def run_STAR_input(wildcards):
    inputs = dict()
    sample_name = wildcards.sample
    fastqs = get_fastqs_for_sample_name(sample_name)
    inputs['fastqs'] = fastqs
    inputs['star_reference'] = os.path.join('references', 'STAR')
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    return inputs


def run_STAR_read_files_in_param(wildcards, input):
    return ' '.join(input.fastqs)


def run_STAR_read_files_command_param(wildcards, input):
    is_gz = list()
    for fastq in input.fastqs:
        is_gz.append(fastq.endswith('.gz'))

    if any(is_gz):
        if not all(is_gz):
            raise Exception('mixed gz and non-gz for sample: {}: {}'
                            .format(wildcards.sample, input.fastqs))

        return '--readFilesCommand zcat'

    return ''


rule run_STAR:
    input:
        unpack(run_STAR_input),
    output:
        bam=os.path.join('bam', '{sample}.bam'),
        idx=os.path.join('bam', '{sample}.bam.bai'),
    log:
        out=os.path.join('bam', 'run_STAR_{sample}_log.out'),
        err=os.path.join('bam', 'run_STAR_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        read_files_in=run_STAR_read_files_in_param,
        read_files_command=run_STAR_read_files_command_param,
        out_prefix=os.path.join('bam', 'STAR_{sample}'),
        star_bam_path=os.path.join('bam', 'STAR_{sample}',
                                   'Aligned.sortedByCoord.out.bam'),
    threads: config['run_STAR_threads']
    resources:
        mem_mb=config['run_STAR_mem_gb'] * 1024,
        time_hours=config['run_STAR_time_hr'],
    shell:
        'mkdir -p {params.out_prefix}'
        ' && {params.conda_wrapper} STAR'
        ' --genomeDir {input.star_reference}'
        ' --sjdbGTFfile {input.gtf}'
        ' --runThreadN {threads}'
        ' --outSAMunmapped Within'
        ' --twopassMode Basic'
        ' --outSAMtype BAM SortedByCoordinate'
        ' --limitBAMsortRAM 24489245788'
        ' --outFilterType BySJout'
        ' --outFilterMultimapNmax 20'
        ' --alignSJoverhangMin 8'
        ' --alignSJDBoverhangMin 1'
        ' --outFilterMismatchNmax 999'
        ' --outFilterMismatchNoverLmax 0.04'
        ' --alignIntronMin 20'
        ' --alignIntronMax 1000000'
        ' --alignMatesGapMax 1000000'
        ' --seedSearchStartLmax 30'
        ' --alignEndsType EndToEnd'
        ' --readFilesIn {params.read_files_in}'
        ' --outFileNamePrefix {params.out_prefix}/'
        ' {params.read_files_command}'
        ' 1> {log.out}'
        ' 2> {log.err}'
        ' && mv {params.star_bam_path} {output.bam}'
        ' 1>> {log.out}'
        ' 2>> {log.err}'
        ' && {params.conda_wrapper} samtools'
        ' index {output.bam}'
        ' 1>> {log.out}'
        ' 2>> {log.err}'


def run_rmats_prep_is_paired_param(wildcards):
    sample_name = wildcards.sample
    fastqs = get_fastqs_for_sample_name(sample_name)
    if len(fastqs) == 2:
        return '--is-paired'
    else:
        return ''


rule run_rmats_prep:
    input:
        bam=os.path.join('bam', '{sample}.bam'),
        gtf=os.path.join('references', config['gtf_name']),
    output:
        dot_rmats=os.path.join('rmats', 'dot_rmats', '{sample}.rmats'),
    log:
        out=os.path.join('rmats', 'run_rmats_prep_{sample}_log.out'),
        err=os.path.join('rmats', 'run_rmats_prep_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'], 'run_rmats_prep_wrapper.py'),
        tmp_dir=os.path.join('rmats', 'tmp_{sample}'),
        out_dir=os.path.join('rmats', 'out_{sample}'),
        b1=os.path.join('rmats', 'b1_{sample}.txt'),
        read_length=config['read_length'],
        is_paired=run_rmats_prep_is_paired_param,
    resources:
        mem_mb=config['run_rmats_prep_mem_gb'] * 1024,
        time_hours=config['run_rmats_prep_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --bam {input.bam}'
        ' --dot-rmats {output.dot_rmats}'
        ' --tmp-dir {params.tmp_dir}'
        ' --out-dir {params.out_dir}'
        ' --b1 {params.b1}'
        ' --gtf {input.gtf}'
        ' --read-length {params.read_length}'
        ' {params.is_paired}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def run_rmats_post_input(wildcards):
    inputs = dict()
    dot_rmats = list()
    for group_key in ['sample_group_1', 'sample_group_2']:
        for name in config[group_key].keys():
            dot_rmats_path = os.path.join(
                'rmats', 'dot_rmats', '{}.rmats'.format(name))
            dot_rmats.append(dot_rmats_path)

    inputs['dot_rmats'] = dot_rmats
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    return inputs


def run_rmats_post_output():
    outputs = dict()
    outputs['mats_jc'] = list()
    outputs['mats_jcec'] = list()
    outputs['fromgtf'] = list()
    outputs['fromgtf_novel_junction'] = list()
    outputs['fromgtf_novel_splice_site'] = list()
    outputs['jc_raw'] = list()
    outputs['jcec_raw'] = list()
    out_dir = os.path.join('rmats', 'out_{cancer_type}')
    for event_type in ALL_EVENT_TYPES:
        mats_jc = '{}.MATS.JC.txt'.format(event_type)
        mats_jcec = '{}.MATS.JCEC.txt'.format(event_type)
        fromgtf = 'fromGTF.{}.txt'.format(event_type)
        fromgtf_novel_junction = 'fromGTF.novelJunction.{}.txt'.format(
            event_type)
        fromgtf_novel_splice_site = 'fromGTF.novelSpliceSite.{}.txt'.format(
            event_type)
        jc_raw = 'JC.raw.input.{}.txt'.format(event_type)
        jcec_raw = 'JCEC.raw.input.{}.txt'.format(event_type)
        outputs['mats_jc'].append(os.path.join(out_dir, mats_jc))
        outputs['mats_jcec'].append(os.path.join(out_dir, mats_jcec))
        outputs['fromgtf'].append(os.path.join(out_dir, fromgtf))
        outputs['fromgtf_novel_junction'].append(os.path.join(
            out_dir, fromgtf_novel_junction))
        outputs['fromgtf_novel_splice_site'].append(os.path.join(
            out_dir, fromgtf_novel_splice_site))
        outputs['jc_raw'].append(os.path.join(out_dir, jc_raw))
        outputs['jcec_raw'].append(os.path.join(out_dir, jcec_raw))

    outputs['b1'] = os.path.join('rmats', 'b1_{cancer_type}.txt')
    outputs['b2'] = os.path.join('rmats', 'b2_{cancer_type}.txt')
    return outputs


def run_rmats_post_b1_bams_param(wildcards):
    bam_paths = list()
    for name in config['sample_group_1']:
        bam_path = os.path.join('bam', '{}.bam'.format(name))
        bam_paths.append(bam_path)

    formatted_paths = ','.join(sorted(bam_paths))
    return '--b1-bams {}'.format(formatted_paths)


def run_rmats_post_b2_bams_param(wildcards):
    if not config.get('sample_group_2'):
        return ''

    bam_paths = list()
    for name in config['sample_group_2']:
        bam_path = os.path.join('bam', '{}.bam'.format(name))
        bam_paths.append(bam_path)

    formatted_paths = ','.join(sorted(bam_paths))
    return '--b2-bams {}'.format(formatted_paths)


def run_rmats_post_dot_rmats_param(wildcards, input):
    formatted_paths = ','.join(input.dot_rmats)
    return '--dot-rmats-files {}'.format(formatted_paths)


rule run_rmats_post:
    input:
        unpack(run_rmats_post_input),
    output:
        **run_rmats_post_output()
    log:
        out=os.path.join('rmats', 'run_rmats_post_{cancer_type}_log.out'),
        err=os.path.join('rmats', 'run_rmats_post_{cancer_type}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'], 'run_rmats_post_wrapper.py'),
        out_dir=os.path.join('rmats', 'out_{cancer_type}'),
        tmp_dir=os.path.join('rmats', 'tmp_{cancer_type}'),
        b1_bams=run_rmats_post_b1_bams_param,
        b2_bams=run_rmats_post_b2_bams_param,
        dot_rmats_files=run_rmats_post_dot_rmats_param,
        read_length=config['read_length'],
    threads: config['run_rmats_post_threads']
    resources:
        mem_mb=config['run_rmats_post_mem_gb'] * 1024,
        time_hours=config['run_rmats_post_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --tmp-dir {params.tmp_dir}'
        ' --out-dir {params.out_dir}'
        ' --b1 {output.b1}'
        ' --b2 {output.b2}'
        ' --gtf {input.gtf}'
        ' --read-length {params.read_length}'
        ' --num-threads {threads}'
        ' {params.dot_rmats_files}'
        ' {params.b1_bams}'
        ' {params.b2_bams}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def generate_psi_matrix_two_group_param(wildcards, input):
    if not config.get('sample_group_2'):
        return ''

    return '--two_group --second_sample_file {}'.format(input.b2)


rule generate_psi_matrix:
    input:
        jc_raw=os.path.join('rmats', 'out_{cancer_type}',
                            'JC.raw.input.{event_type}.txt'),
        from_gtf=os.path.join('rmats', 'out_{cancer_type}',
                              'fromGTF.{event_type}.txt'),
        from_gtf_novel_junction=os.path.join(
            'rmats', 'out_{cancer_type}',
            'fromGTF.novelJunction.{event_type}.txt'),
        from_gtf_novel_splice=os.path.join(
            'rmats', 'out_{cancer_type}',
            'fromGTF.novelSpliceSite.{event_type}.txt'),
        b1=os.path.join('rmats', 'b1_{cancer_type}.txt'),
        b2=os.path.join('rmats', 'b2_{cancer_type}.txt'),
    output:
        merged=os.path.join('rMATS_matrix',
                            'merged_matrix.{event_type}.{cancer_type}.txt'),
        annotation=os.path.join('rMATS_matrix',
                                'AS_annotation.{event_type}.{cancer_type}.txt'),
    log:
        out=os.path.join(
            'rMATS_matrix',
            'generate_psi_matrix_{cancer_type}_{event_type}_log.out'),
        err=os.path.join(
            'rMATS_matrix',
            'generate_psi_matrix_{cancer_type}_{event_type}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'], 'generate_PSI_matrix.py'),
        rmats_out_dir=os.path.join('rmats', 'out_{cancer_type}'),
        out_dir='rMATS_matrix',
        two_group=generate_psi_matrix_two_group_param,
    resources:
        mem_mb=config['generate_psi_matrix_mem_gb'] * 1024,
        time_hours=config['generate_psi_matrix_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' {params.rmats_out_dir}'
        ' {params.out_dir}'
        ' --na_filter'
        ' --coverage_filter'
        ' --event {wildcards.event_type}'
        ' --tissue {wildcards.cancer_type}'
        ' --event_annotation'
        ' --sample {input.b1}'
        ' {params.two_group}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule filter_psi_matrix:
    input:
        merged=os.path.join('rMATS_matrix',
                            'merged_matrix.{event_type}.{cancer_type}.txt'),
    output:
        filtered=os.path.join('rMATS_matrix',
                              'filtered_psi.matrix.{cancer_type}.{event_type}.csv'),
    log:
        out=os.path.join(
            'rMATS_matrix',
            'filter_psi_matrix_{cancer_type}_{event_type}_log.out'),
        err=os.path.join(
            'rMATS_matrix',
            'filter_psi_matrix_{cancer_type}_{event_type}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'], 'filter_PSI_matrix.py'),
        out_dir='rMATS_matrix',
    resources:
        mem_mb=config['filter_psi_matrix_mem_gb'] * 1024,
        time_hours=config['filter_psi_matrix_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' {input.merged}'
        ' {params.out_dir}'
        ' --cancer_type {wildcards.cancer_type}'
        ' --rna'
        ' 1> {log.out}'
        ' 2> {log.err}'


def generate_psi_matrix_from_mats_input(wildcards):
    inputs = dict()
    rmats_dir = os.path.join('rmats', 'out_{}'.format(wildcards.cancer_type))
    for event_type in USED_EVENT_TYPES:
        inputs[event_type] = os.path.join(rmats_dir,
                                          '{}.MATS.JC.txt'.format(event_type))

    return inputs


def generate_psi_matrix_from_mats_output():
    outputs = dict()
    out_dir = 'rMATS_matrix'
    for event_type in USED_EVENT_TYPES:
        outputs[event_type] = os.path.join(
            out_dir, '{{cancer_type}}.{}.MATS_df.txt'.format(event_type))

    return outputs


rule generate_psi_matrix_from_mats:
    input:
        unpack(generate_psi_matrix_from_mats_input),
    output:
        **generate_psi_matrix_from_mats_output()
    log:
        out=os.path.join('rMATS_matrix',
                         'generate_psi_matrix_from_mats_{cancer_type}_log.out'),
        err=os.path.join('rMATS_matrix',
                         'generate_psi_matrix_from_mats_{cancer_type}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'], 'generate_PSI_matrix_from_MATS.py'),
        rmats_dir = os.path.join('rmats', 'out_{cancer_type}'),
        matrix_dir = 'rMATS_matrix',
    resources:
        mem_mb=config['generate_psi_matrix_from_mats_mem_gb'] * 1024,
        time_hours=config['generate_psi_matrix_from_mats_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' {params.rmats_dir}'
        ' {params.matrix_dir}'
        ' {wildcards.cancer_type}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def quantify_expression_single_end_param(wildcards):
    if config.get('sample_group_2'):
        return ''

    read_length = config['read_length']
    std_dev_of_length = '10'
    return '--single -l {} -s {}'.format(read_length, std_dev_of_length)


def quantify_expression_r1_view_command_param(wildcards, input):
    pipe_command = '| samtools bam2fq -'
    if config.get('sample_group_2'):
        view_command = 'samtools view -uf64 {}'.format(input.bam)
        return '{} {}'.format(view_command, pipe_command)

    view_command = 'samtools view -F 0x4 {} -b'.format(input.bam)
    return '{} {}'.format(view_command, pipe_command)


def quantify_expression_r2_view_command_param(wildcards, input):
    # make an empty tmp_r2 when there is only 1 sample group
    if not config.get('sample_group_2'):
        return 'echo ""'

    pipe_command = '| samtools bam2fq -'
    view_command = 'samtools view -uf128 {}'.format(input.bam)
    return '{} {}'.format(view_command, pipe_command)


def quantify_expression_kallisto_r2_param(wildcards, output):
    if not config.get('sample_group_2'):
        return ''

    return output.tmp_r2


rule quantify_expression:
    input:
        bam=os.path.join('bam', '{sample}.bam'),
        idx=os.path.join('references', 'homo_sapiens', 'transcriptome.idx'),
    output:
        abundance=os.path.join('expression', '{sample}_quant', 'abundance.tsv'),
        tmp_r1=temp(os.path.join('expression', '{sample}_quant', '{sample}_R1')),
        tmp_r2=temp(os.path.join('expression', '{sample}_quant', '{sample}_R2')),
    log:
        out=os.path.join('expression', 'quantify_expression_{sample}_log.out'),
        err=os.path.join('expression', 'quantify_expression_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        out_dir=os.path.join('expression', '{sample}_quant'),
        read_length=config['read_length'],
        single_end=quantify_expression_single_end_param,
        r1_view_command=quantify_expression_r1_view_command_param,
        r2_view_command=quantify_expression_r2_view_command_param,
        kallisto_r2=quantify_expression_kallisto_r2_param,
    resources:
        mem_mb=config['quantify_expression_mem_gb'] * 1024,
        time_hours=config['quantify_expression_time_hr'],
    shell:
        '{params.conda_wrapper}'
        ' {params.r1_view_command}'
        ' > {output.tmp_r1}'
        ' 2> {log.err}'
        ' && {params.conda_wrapper}'
        ' {params.r2_view_command}'
        ' > {output.tmp_r2}'
        ' 2>> {log.err}'
        ' && {params.conda_wrapper}'
        ' kallisto quant'
        ' -i {input.idx}'
        ' -o {params.out_dir}'
        ' {output.tmp_r1}'
        ' {params.kallisto_r2}'
        ' {params.single_end}'
        ' 1> {log.out}'
        ' 2>> {log.err}'


def convert_and_merge_kallisto_run_input(wildcards):
    inputs = dict()
    inputs['to_genes'] = os.path.join('references', 'homo_sapiens',
                                      'transcripts_to_genes.txt'),
    for group_key in ['sample_group_1', 'sample_group_2']:
        for name in config[group_key].keys():
            inputs[name] = os.path.join(
                'expression', '{}_quant'.format(name), 'abundance.tsv')

    return inputs


rule convert_and_merge_kallisto_run:
    input:
        unpack(convert_and_merge_kallisto_run_input),
    output:
        merged=os.path.join('expression',
                            '{cancer_type}.merged.exp.kallisto.txt'),
    log:
        out=os.path.join('expression',
                         'convert_and_merge_kallisto_run_{cancer_type}_log.out'),
        err=os.path.join('expression',
                         'convert_and_merge_kallisto_run_{cancer_type}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join(config['scripts_dir'],
                            'convert_and_merge_kallisto_run.py'),
        in_dir='expression',
        out_dir='expression',
    resources:
        mem_mb=config['convert_and_merge_kallisto_run_mem_gb'] * 1024,
        time_hours=config['convert_and_merge_kallisto_run_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' {params.in_dir}'
        ' {params.out_dir}'
        ' {wildcards.cancer_type}'
        ' --transcripts-to-genes {input.to_genes}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def generate_SPARKS_object_input(wildcards):
    inputs = dict()
    two_group = bool(config.get('sample_group_2'))
    inputs['exp_matrix'] = os.path.join(
        'expression', '{}.merged.exp.kallisto.txt'.format(wildcards.cancer_type))
    for event_type in USED_EVENT_TYPES:
        inputs['psi_{}'.format(event_type)] = os.path.join(
            'rMATS_matrix',
            'filtered_psi.matrix.{}.{}.csv'.format(wildcards.cancer_type,
                                                   event_type))
        inputs['anno_{}'.format(event_type)] = os.path.join(
            'rMATS_matrix',
            'AS_annotation.{}.{}.txt'.format(event_type,
                                             wildcards.cancer_type))
        if two_group:
            inputs['mats_{}'.format(event_type)] = os.path.join(
                'rMATS_matrix', '{}.{}.MATS_df.txt'.format(wildcards.cancer_type,
                                                           event_type))

    return inputs


def generate_SPARKS_object_psi_param(wildcards, input):
    args = list()
    for event_type in USED_EVENT_TYPES:
        args.append('--psi_{}'.format(event_type))
        args.append(input['psi_{}'.format(event_type)])

    return ' '.join(args)


def generate_SPARKS_object_anno_param(wildcards, input):
    args = list()
    for event_type in USED_EVENT_TYPES:
        args.append('--exon_anno_{}'.format(event_type))
        args.append(input['anno_{}'.format(event_type)])

    return ' '.join(args)


def generate_SPARKS_object_mats_param(wildcards, input):
    if not config.get('sample_group_2'):
        return ''

    args = list()
    for event_type in USED_EVENT_TYPES:
        args.append('--MATS_{}'.format(event_type))
        args.append(input['mats_{}'.format(event_type)])

    return ' '.join(args)


def generate_SPARKS_object_sparks_lib_param(wildcards):
    if not config['analysis_mode']:
        return ''

    # TODO
    lib_path = os.path.join(
        '/mnt', 'isilon', 'xing_lab', 'aspera', 'yangt3', 'ENCODE', 'SPARKS',
        'ENCODE_library.KD_and_KO.count_20.rds')
    return '--SPARKS_library {}'.format(lib_path)


rule generate_SPARKS_object:
    input:
        unpack(generate_SPARKS_object_input),
    output:
        rds='{cancer_type}.SPARKS.rds',
    log:
        out='generate_SPARKS_object_{cancer_type}_log.out',
        err='generate_SPARKS_object_{cancer_type}_log.err',
    params:
        conda_wrapper=config['conda_wrapper'],
        script='generate_SPARKS_object.R',
        sparks_lib=generate_SPARKS_object_sparks_lib_param,
        psi=generate_SPARKS_object_psi_param,
        anno=generate_SPARKS_object_anno_param,
        mats=generate_SPARKS_object_mats_param,
    resources:
        mem_mb=config['generate_SPARKS_object_mem_gb'] * 1024,
        time_hours=config['generate_SPARKS_object_time_hr'],
    shell:
        '{params.conda_wrapper} Rscript'
        ' {params.script}'
        ' --exp {input.exp_matrix}'
        ' --cancer_type {wildcards.cancer_type}'
        ' {params.sparks_lib}'
        ' {params.psi}'
        ' {params.anno}'
        ' {params.mats}'
        ' 1> {log.out}'
        ' 2> {log.err}'
