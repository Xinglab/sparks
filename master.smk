from glob import glob
import os
configfile: "./config.yaml"

# path to references
REPO_DIR = "~/xinglab/tools/SPARKS/"
REF_FASTA = "~/xinglab/references/grch37/ucsc.hg19.fasta"
KAL_INDEX_FILE = "/home/yangt3/xinglab/references/kallisto/homo_sapiens/transcriptome.idx"
GENCODE_GTF = "/home/yangt3/aspera/GTEx/references/37.chr/gencode.v26lift37.annotation.gtf"
SPARKS_LIB_FILE = "/home/yangt3/xinglab/SPARKS/library/Clean_minfiltered.combined.library.rds"
WD = os.path.dirname(os.path.realpath(__file__))

# event type definition
event_types = ["SE", "A3SS", "A5SS", "RI"]
lm_result_type = ["beta", "pval"]

# configuration
# define study
study = config['study']
print(study)

# define read length
READ_LENGTH = config['read_length']
print("Read Length : %s" % READ_LENGTH)

# rule order definition based on config

# define samples - explicitly declared from config rather than implicit glob for all
# and configure (case vs. control) vs. run w/o statistical part

# extract identifiers
case_identifiers = config['sample_definition']['case']
ctrl_identifiers = config['sample_definition']['control']


# configure paired end vs. single end reads by parsing sample names 
if config['paired_end']:  # if paired end 
    ruleorder: run_STAR_paired_end > run_STAR_single_end

    ruleorder: quantify_expression_paired_end > quantify_expression_single_end

    # define char for rmats run
    READ_PAIRED = 'paired'

    print("Paired end mode selected for reads")
    # samples, = glob_wildcards("fastq/{sample}_1.fastq.gz")
    # identify sample names
    case_files = [filename
                  for filenames in [glob('fastq/%s.*_1.fastq.gz' % x) for x in case_identifiers]
                  for filename in filenames]
    ctrl_files = [filename
                  for filenames in [glob('fastq/%s.*_1.fastq.gz' % x) for x in ctrl_identifiers]
                  for filename in filenames]

    # parse sample names
    case_samples = [x.split('/')[-1].split('_1.fastq')[0] for x in case_files]
    ctrl_samples = [x.split('/')[-1].split('_1.fastq')[0] for x in ctrl_files]

else:  # if single end 
    ruleorder: run_STAR_single_end > run_STAR_paired_end

    ruleorder: quantify_expression_single_end > quantify_expression_paired_end

    # define char for rmats run
    READ_PAIRED = 'single'

    print("Single end mode selected for reads")
    # samples, = glob_wildcards("fastq/{sample}.fastq.gz")
    # identify sample names
    case_files = [filename
                  for filenames in [glob('fastq/%s.*.fastq.gz' % x) for x in case_identifiers]
                  for filename in filenames]
    ctrl_files = [filename
                  for filenames in [glob('fastq/%s.*.fastq.gz' % x) for x in ctrl_identifiers]
                  for filename in filenames]

    # parse sample names
    case_samples = [x.split('/')[-1].split('.fastq')[0] for x in case_files]
    ctrl_samples = [x.split('/')[-1].split('.fastq')[0] for x in ctrl_files]

samples = case_samples + ctrl_samples
print(samples)

# configure rMATS stats mode 
if config['stat_mode']:
    ruleorder: gather_case_ctrl_samples_for_stat > gather_case_ctrl_samples_wo_stat
else:
    ruleorder: gather_case_ctrl_samples_wo_stat > gather_case_ctrl_samples_for_stat


# configure clip
if config['clip_mode']:
    print("CLIP analysis mode enabled")
    if config['stat_mode']:
        if config['analysis_mode']:
            ruleorder: generate_SPARKS_full_with_analysis > generate_SPARKS_full > generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_full_without_stat > generate_SPARKS_naive
        else: 
            ruleorder: generate_SPARKS_full > generate_SPARKS_full_with_analysis > generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_full_without_stat > generate_SPARKS_naive
    else:  # no stat mode
        ruleorder: generate_SPARKS_full_without_stat > generate_SPARKS_full > generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_naive
else:
    print("CLIP analysis mode disabled")
    if config['stat_mode']:
        if config['analysis_mode']:
            ruleorder: generate_SPARKS_full_with_stat_without_clip_with_analysis > generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_naive > generate_SPARKS_full >  generate_SPARKS_full_without_stat
        else:
            ruleorder: generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_full_with_stat_without_clip_with_analysis > generate_SPARKS_naive > generate_SPARKS_full >  generate_SPARKS_full_without_stat            
    else:  # no stat mode
        ruleorder: generate_SPARKS_naive > generate_SPARKS_full_with_stat_without_clip > generate_SPARKS_full >  generate_SPARKS_full_without_stat


# make directory for lgos 
if not os.path.exists('./logs/'):
    os.makedirs('./logs/')

##### RULES ##### 
rule all:
    input:
        # "{sample}.TumorSplicingExperiment.rds"
        "{cancer_type}.done"

# Run alignment
rule run_STAR_paired_end:
    input:
        expand("fastq/{{sample}}_{read}.fastq.gz",read=['1', '2'])
    output:
        "bam/{sample}.bam"
    params:
        # cluster = "-l h_vmem=64G"
        cluster="--mem=64G -t 12:00:00 --cpus-per-task=8"
    shell:
        """
        # module load slurm
        # module load STAR
        export LD_LIBRARY_PATH=/cm/shared/easybuild/software/GCCcore/10.3.0/lib64/:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64:${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}
        module load STAR 
        # module load SAMtools
        mkdir --parents bam/STAR_{wildcards.sample}/ 
        STAR \
        --genomeDir /home/yangt3/xinglab/references/STAR \
        --sjdbGTFfile {GENCODE_GTF} \
        --runThreadN 8 \
        --outSAMunmapped Within \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 24489245788 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --seedSearchStartLmax 30 \
        --alignEndsType EndToEnd \
        --readFilesIn fastq/{wildcards.sample}_1.fastq.gz fastq/{wildcards.sample}_2.fastq.gz \
        --outFileNamePrefix bam/STAR_{wildcards.sample}/ \
        --readFilesCommand zcat
        mv bam/STAR_{wildcards.sample}/Aligned.sortedByCoord.out.bam bam/{wildcards.sample}.bam
        /home/yangt3/xinglab/tools/samtools-1.12/samtools index bam/{wildcards.sample}.bam
        """

rule run_STAR_single_end:
    input:
        expand("fastq/{{sample}}.fastq.gz")
    output:
        "bam/{sample}.bam"
    params:
        # cluster = "-l h_vmem=64G"
        cluster="--mem=64G -t 12:00:00 --cpus-per-task=8"
    shell:
        """
        # module load slurm
        # module load STAR
        export LD_LIBRARY_PATH=/cm/shared/easybuild/software/GCCcore/10.3.0/lib64/:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64:/cm/shared/easybuild/software/GCCcore/10.3.0/lib64/:${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}
        module load STAR 
        # module load SAMtools
        mkdir --parents bam/STAR_{wildcards.sample}/ 
        STAR \
        --genomeDir /home/yangt3/xinglab/references/STAR \
        --sjdbGTFfile {GENCODE_GTF} \
        --runThreadN 8 \
        --outSAMunmapped Within \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 24489245788 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --seedSearchStartLmax 30 \
        --alignEndsType EndToEnd \
        --readFilesIn fastq/{wildcards.sample}.fastq.gz \
        --outFileNamePrefix bam/STAR_{wildcards.sample}/ \
        --readFilesCommand zcat
        mv bam/STAR_{wildcards.sample}/Aligned.sortedByCoord.out.bam bam/{wildcards.sample}.bam
        /home/yangt3/xinglab/tools/samtools-1.12/samtools index bam/{wildcards.sample}.bam
        """

# generate sample list for rMATS input
rule gather_case_ctrl_samples_for_stat:
    input:
        case_bams = expand('bam/{sample}.bam',sample=case_samples),
        ctrl_bams = expand('bam/{sample}.bam',sample=ctrl_samples)
    output:
        case_sample_list = "rMATS_stats/sample_case.txt",
        ctrl_sample_list = "rMATS_stats/sample_ctrl.txt"
    params:
        cluster = "--mem=1G -t 01:00:00"

    run:
        import os
        if not os.path.exists('rMATS_stats'):
            os.mkdir('rMATS_stats')
        cwd = os.getcwd()
        with open(output.case_sample_list, 'w') as f:
            output_sample_list = [ "%s/%s" % (cwd, x) for x in input.case_bams ]
            output_string = ",".join(sorted(output_sample_list))
            f.write(output_string)
        with open(output.ctrl_sample_list, 'w') as f:
            output_sample_list = [ "%s/%s" % (cwd, x) for x in input.ctrl_bams ]
            output_string = ",".join(sorted(output_sample_list))
            f.write(output_string)

rule gather_case_ctrl_samples_wo_stat:
    input:
        sample_bams = expand('bam/{sample}.bam',sample=samples),
    output:
        sample_list = "rMATS_naive/sample_list.txt"
    params:
        cluster = "--mem=1G -t 01:00:00"

    run:
        import os
        if not os.path.exists('rMATS_naive'):
            os.mkdir('rMATS_naive')
        cwd = os.getcwd()
        with open(output.sample_list, 'w') as f:
            output_sample_list = [ "%s/%s" % (cwd, x) for x in input.sample_bams ]
            output_string = ",".join(output_sample_list)
            f.write(output_string)

# Run rMATS
rule run_rmats_with_stat:
    input:
        expand('bam/{sample}.bam',sample=samples),
        "rMATS_stats/sample_case.txt",
        "rMATS_stats/sample_ctrl.txt"

    output:
        expand("rMATS_stats/JC.raw.input.{event_type}.txt",event_type=event_types),
        expand("rMATS_stats/{event_type}.MATS.JC.txt",event_type=event_types),
        expand("rMATS_stats/fromGTF.{event_type}.txt",event_type=event_types)

    params:
        # cluster = "-l h_vmem=64G"
        cluster="--mem=64G -t 24:00:00 --cpus-per-task=6"
    shell:
        """
        export LD_LIBRARY_PATH=/home/yangt3/xinglab/tools/gsl-2.6/lib:/usr/lib64/:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64/:/cm/shared/apps/openblas/0.3.7/lib/:/home/yangt3/xinglab/tools/OpenBLAS/:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64/:/home/yangt3/xinglab/tools/rmats-turbo/conda_envs/rmats/lib/
        module load lapack/3.9.0
        
        /home/yangt3/xinglab/tools/rmats-turbo-4.1.1/run_rmats \
        --b1 $PWD/rMATS_stats/sample_ctrl.txt \
        --b2 $PWD/rMATS_stats/sample_case.txt \
        --gtf {GENCODE_GTF} \
        --od $PWD/rMATS_stats/ \
        -t {READ_PAIRED} \
        --readLength {READ_LENGTH} \
        --variable-read-length \
        --nthread 6 \
        --tmp $PWD/rMATS_stats/tmp \
        --novelSS 
        """

rule run_rmats_without_stat:
    input:
        expand('bam/{sample}.bam',sample=samples),
        "rMATS_naive/sample_list.txt",

    output:
        expand("rMATS_naive/JC.raw.input.{event_type}.txt",event_type=event_types),
        expand("rMATS_naive/fromGTF.{event_type}.txt",event_type=event_types)

    params:
        # cluster = "-l h_vmem=64G"
        cluster="--mem=64G -t 24:00:00 --cpus-per-task=6"
    shell:
        """
        export LD_LIBRARY_PATH=/home/yangt3/xinglab/tools/gsl-2.6/lib:/usr/lib64/:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64/:/cm/shared/apps/openblas/0.3.7/lib/:/home/yangt3/xinglab/tools/OpenBLAS/:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64/:/home/yangt3/xinglab/tools/rmats-turbo/conda_envs/rmats/lib/
        module load lapack/3.9.0
        
        /home/yangt3/xinglab/tools/rmats-turbo/run_rmats \
        --b1 $PWD/rMATS_naive/sample_list.txt \
        --gtf {GENCODE_GTF} \
        --od $PWD/rMATS_naive/ \
        -t {READ_PAIRED} \
        --readLength {READ_LENGTH} \
        --variable-read-length \
        --nthread 6 \
        --tmp $PWD/rMATS_naive/tmp \
        --novelSS \
        --statoff
        """

# Generate PSI matrix 
rule psi_matrix_generation_without_stat:
    input:
        expand("rMATS_naive/JC.raw.input.{event_type}.txt",event_type=event_types)
    output:
        "rMATS_naive_matrix/merged_matrix.{event_type}.{cancer_type}.txt",
        "rMATS_naive_matrix/AS_annotation.{event_type}.{cancer_type}.txt"
    params:
        # cluster = "-l h_vmem=8G"
        cluster="--mem=8G -t 06:00:00"
    shell:
        """
        # python /home/yangt3/xinglab/gtex_sqtl/gtex_sqtl/scripts/summary/generate_PSI_matrix.py \

        python {REPO_DIR}/aux_scripts/generate_PSI_matrix.py \
        $PWD/rMATS_naive/ \
        $PWD/rMATS_naive_matrix/ \
        --na_filter \
        --coverage_filter \
        --event {wildcards.event_type} \
        --tissue {wildcards.cancer_type} \
        --event_annotation \
        --sample $PWD/rMATS_naive/sample_list.txt
        """

rule psi_matrix_generation_with_stat:
    input:
        expand("rMATS_stats/JC.raw.input.{event_type}.txt",event_type=event_types)
    output:
        "rMATS_stats_matrix/merged_matrix.{event_type}.{cancer_type}.txt",
        "rMATS_stats_matrix/AS_annotation.{event_type}.{cancer_type}.txt"
    params:
        # cluster = "-l h_vmem=8G"
        cluster="--mem=8G -t 06:00:00"
    shell:
        """
        # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64

        python {REPO_DIR}/aux_scripts/generate_PSI_matrix.py \
        $PWD/rMATS_stats/ \
        $PWD/rMATS_stats_matrix/ \
        --na_filter \
        --coverage_filter \
        --event {wildcards.event_type} \
        --tissue {wildcards.cancer_type} \
        --event_annotation \
        --two_group \
        --second_sample_file $PWD/rMATS_stats/sample_case.txt \
        --sample $PWD/rMATS_stats/sample_ctrl.txt
        """


rule filter_psi_matrix_without_stat:
    input:
        psi_matrix="rMATS_naive_matrix/merged_matrix.{event_type}.{cancer_type}.txt"
    output:
        "rMATS_naive_matrix/filtered_psi.matrix.{cancer_type}.{event_type}.csv"
    params:
        # cluster = "-l h_vmem=32G"
        cluster="--mem=32G -t 10:00:00"
    shell:
        """
        python /{REPO_DIR}/aux_scripts/filter_PSI_matrix.py {input.psi_matrix} \
        $PWD/rMATS_naive_matrix/ --cancer_type {wildcards.cancer_type} --rna
        """

rule filter_psi_matrix_with_stat:
    input:
        psi_matrix="rMATS_stats_matrix/merged_matrix.{event_type}.{cancer_type}.txt"
    output:
        "rMATS_stats_matrix/filtered_psi.matrix.{cancer_type}.{event_type}.csv"
    params:
        # cluster = "-l h_vmem=32G"
        cluster="--mem=32G -t 10:00:00"
    shell:
        """
        python {REPO_DIR}/aux_scripts/filter_PSI_matrix.py {input.psi_matrix} \
        $PWD/rMATS_stats_matrix/ --cancer_type {wildcards.cancer_type} --rna
        """

# Gather MATS result for SPARKS object input
rule gather_MATS_result:
    input:
        expand("rMATS_stats/{event_type}.MATS.JC.txt",event_type=event_types)
    output:
        expand("rMATS_stats_matrix/{{cancer_type}}.{event_type}.MATS_df.txt",event_type=event_types)
    params:
        # cluster = "-l h_vmem=8G"
        cluster="--mem=8G -t 06:00:00"
    shell:
        """
        python {REPO_DIR}/aux_scripts/generate_PSI_matrix_from_MATS.py \
        $PWD/rMATS_stats/ \
        $PWD/rMATS_stats_matrix/ \
        {wildcards.cancer_type}

        """

# Quantify expression

rule quantify_expression_single_end:
    input:
        'bam/{sample}.bam'
    output:
        'expression/{sample}_quant/abundance.tsv',
        temp("expression/{sample}_quant/{sample}_R1")
    params:
        # cluster = "-l h_vmem=8G"
        cluster="--mem=8G -t 06:00:00"
    run:
        shell("mkdir -p expression/{wildcards.sample}_quant")
        # shell("mkfifo expression/{wildcards.sample}_quant/{wildcards.sample}_R1 && /home/yangt3/xinglab/tools/samtools-1.12/samtools view -F 0x4 {input} | /home/yangt3/xinglab/tools/samtools-1.12/samtools bam2fq - > expression/{wildcards.sample}_quant/{wildcards.sample}_R1 & \
            # kallisto quant -i {KAL_INDEX_FILE} -o expression/{wildcards.sample}_quant expression/{wildcards.sample}_quant/{wildcards.sample}_R1 --single -l {READ_LENGTH} -s 15")
        shell("/home/yangt3/xinglab/tools/samtools-1.12/samtools view -F 0x4 {input} -b | /home/yangt3/xinglab/tools/samtools-1.12/samtools bam2fq - > expression/{wildcards.sample}_quant/{wildcards.sample}_R1 & \
            kallisto quant -i {KAL_INDEX_FILE} -o expression/{wildcards.sample}_quant expression/{wildcards.sample}_quant/{wildcards.sample}_R1 --single -l {READ_LENGTH} -s 10")

rule quantify_expression_paired_end:
    input:
        'bam/{sample}.bam'
    output:
        'expression/{sample}_quant/abundance.tsv',
        temp("expression/{sample}_quant/{sample}_R1"),
        temp("expression/{sample}_quant/{sample}_R2")
    params:
        # cluster = "-l h_vmem=8G"
        cluster="--mem=8G -t 06:00:00"
    run:
        shell("mkdir -p expression/{wildcards.sample}_quant")
        shell("mkfifo expression/{wildcards.sample}_quant/{wildcards.sample}_R1 && /home/yangt3/xinglab/tools/samtools-1.12/samtools view -uf64 {input} | /home/yangt3/xinglab/tools/samtools-1.12/samtools bam2fq - > expression/{wildcards.sample}_quant/{wildcards.sample}_R1 & \
            mkfifo expression/{wildcards.sample}_quant/{wildcards.sample}_R2 && /home/yangt3/xinglab/tools/samtools-1.12/samtools view -uf128 {input} | /home/yangt3/xinglab/tools/samtools-1.12/samtools bam2fq - > expression/{wildcards.sample}_quant/{wildcards.sample}_R2 & \
            kallisto quant -i {KAL_INDEX_FILE} -o expression/{wildcards.sample}_quant expression/{wildcards.sample}_quant/{wildcards.sample}_R1 expression/{wildcards.sample}_quant/{wildcards.sample}_R2")

# merge the expression files into matrix for downstream analysis 
rule merge_exp_into_matrix:
    input:
        expand('expression/{sample}_quant/abundance.tsv',sample=samples)
    output:
        "{cancer_type}.merged.exp.kallisto.txt"
    params:
        # cluster = "-l h_vmem=32G"
        cluster="--mem=32G -t 06:00:00"
    shell:
        "python {REPO_DIR}/aux_scripts/convert_and_merge_kallisto_run.py $PWD/expression/ $PWD {wildcards.cancer_type}"

# run CLIP analysis by intersecting splicing regions and CLIP peak bed files 
rule intersect_CLIP:
    input:
        # expand("hotspot_analysis_{event_type}/linear_model/lm_{result_type}.{batch_id}.csv", result_type = lm_result_type, allow_missing=True)
        "rMATS_stats_matrix/filtered_psi.matrix.{cancer_type}.{event_type}.csv"
    output:
        "CLIP/{cancer_type}.CLIP_overlap.{event_type}.tsv"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        # module load BEDTools
        export PATH=$PATH:/home/yangt3/xinglab/tools/bedtools2/bin/
        python {REPO_DIR}/aux_scripts/find_CLIP_peak_overlap.py {input} $PWD/CLIP/ {wildcards.cancer_type} {wildcards.event_type}
        """

# Generate SPARKS objects 
rule generate_SPARKS_full_with_analysis:
    input:
        psi_df_list=expand("rMATS_stats_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt",
        clip_result=expand("CLIP/{{cancer_type}}.CLIP_overlap.{event_type}.tsv",event_type=event_types),
        mats_result=expand("rMATS_stats_matrix/{{cancer_type}}.{event_type}.MATS_df.txt",event_type=event_types)
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00 --cpus-per-task=3"
    shell:
        """
        export LD_LIBRARY_PATH=/cm/shared/easybuild/software/GCCcore/10.3.0/lib64:/cm/shared/easybuild/software/bzip2/1.0.8-GCCcore-10.3.0/lib64:/cm/local/apps/gcc/9.2.0/lib64:${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64:/cm/shared/apps_chop/readline/8.1/lib
        module load R
        # module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --CLIP_SE CLIP/{wildcards.cancer_type}.CLIP_overlap.SE.tsv \
        --CLIP_A3SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A3SS.tsv \
        --CLIP_A5SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A5SS.tsv \
        --CLIP_RI CLIP/{wildcards.cancer_type}.CLIP_overlap.RI.tsv \
        --exon_anno_SE rMATS_stats_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rMATS_stats_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rMATS_stats_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_stats_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt \
        --MATS_SE rMATS_stats_matrix/{wildcards.cancer_type}.SE.MATS_df.txt \
        --MATS_A3SS rMATS_stats_matrix/{wildcards.cancer_type}.A3SS.MATS_df.txt \
        --MATS_A5SS rMATS_stats_matrix/{wildcards.cancer_type}.A5SS.MATS_df.txt \
        --MATS_RI rMATS_stats_matrix/{wildcards.cancer_type}.RI.MATS_df.txt \
        --SPARKS_library {SPARKS_LIB_FILE}

        touch {wildcards.cancer_type}.done

        """


rule generate_SPARKS_full:
    input:
        psi_df_list=expand("rMATS_stats_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt",
        clip_result=expand("CLIP/{{cancer_type}}.CLIP_overlap.{event_type}.tsv",event_type=event_types),
        mats_result=expand("rMATS_stats_matrix/{{cancer_type}}.{event_type}.MATS_df.txt",event_type=event_types)
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        export LD_LIBRARY_PATH=/cm/shared/easybuild/software/GCCcore/10.3.0/lib64:/cm/shared/easybuild/software/bzip2/1.0.8-GCCcore-10.3.0/lib64:/cm/local/apps/gcc/9.2.0/lib64:${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64:/cm/shared/apps_chop/readline/8.1/lib
        module load R
        # module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --CLIP_SE CLIP/{wildcards.cancer_type}.CLIP_overlap.SE.tsv \
        --CLIP_A3SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A3SS.tsv \
        --CLIP_A5SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A5SS.tsv \
        --CLIP_RI CLIP/{wildcards.cancer_type}.CLIP_overlap.RI.tsv \
        --exon_anno_SE rMATS_stats_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rMATS_stats_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rMATS_stats_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_stats_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt \
        --MATS_SE rMATS_stats_matrix/{wildcards.cancer_type}.SE.MATS_df.txt \
        --MATS_A3SS rMATS_stats_matrix/{wildcards.cancer_type}.A3SS.MATS_df.txt \
        --MATS_A5SS rMATS_stats_matrix/{wildcards.cancer_type}.A5SS.MATS_df.txt \
        --MATS_RI rMATS_stats_matrix/{wildcards.cancer_type}.RI.MATS_df.txt 

        touch {wildcards.cancer_type}.done

        """


rule generate_SPARKS_full_with_stat_without_clip_with_analysis:
    input:
        psi_df_list=expand("rmats_stats_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt",
        mats_result=expand("rmats_stats_matrix/{{cancer_type}}.{event_type}.MATS_df.txt",event_type=event_types),
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64
        module load R
        module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --exon_anno_SE rmats_stats_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rmats_stats_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rmats_stats_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_stats_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt \
        --MATS_SE rmats_stats_matrix/{wildcards.cancer_type}.SE.MATS_df.txt \
        --MATS_A3SS rmats_stats_matrix/{wildcards.cancer_type}.A3SS.MATS_df.txt \
        --MATS_A5SS rmats_stats_matrix/{wildcards.cancer_type}.A5SS.MATS_df.txt \
        --MATS_RI rMATS_stats_matrix/{wildcards.cancer_type}.RI.MATS_df.txt \
        --SPARKS_library {SPARKS_LIB_FILE}
    
        touch {wildcards.cancer_type}.done

        """

rule generate_SPARKS_full_with_stat_without_clip:
    input:
        psi_df_list=expand("rmats_stats_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt",
        mats_result=expand("rmats_stats_matrix/{{cancer_type}}.{event_type}.MATS_df.txt",event_type=event_types),
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64
        module load R
        module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rmats_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_stats_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --exon_anno_SE rmats_stats_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rmats_stats_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rmats_stats_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_stats_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt \
        --MATS_SE rmats_stats_matrix/{wildcards.cancer_type}.SE.MATS_df.txt \
        --MATS_A3SS rmats_stats_matrix/{wildcards.cancer_type}.A3SS.MATS_df.txt \
        --MATS_A5SS rmats_stats_matrix/{wildcards.cancer_type}.A5SS.MATS_df.txt \
        --MATS_RI rMATS_stats_matrix/{wildcards.cancer_type}.RI.MATS_df.txt 
    
        touch {wildcards.cancer_type}.done

        """

rule generate_SPARKS_full_without_stat:
    input:
        psi_df_list=expand("rMATS_naive_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt",
        clip_result=expand("CLIP/{{cancer_type}}.CLIP_overlap.{event_type}.tsv",event_type=event_types)
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64
        module load R
        module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --CLIP_SE CLIP/{wildcards.cancer_type}.CLIP_overlap.SE.tsv \
        --CLIP_A3SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A3SS.tsv \
        --CLIP_A5SS CLIP/{wildcards.cancer_type}.CLIP_overlap.A5SS.tsv \
        --CLIP_RI CLIP/{wildcards.cancer_type}.CLIP_overlap.RI.tsv \
        --exon_anno_SE rMATS_naive_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rMATS_naive_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rMATS_naive_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_naive_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt

        touch {wildcards.cancer_type}.done

        """

rule generate_SPARKS_naive:
    input:
        psi_df_list=expand("rMATS_naive_matrix/filtered_psi.matrix.{{cancer_type}}.{event_type}.csv",event_type=event_types),
        exp_matrix="{cancer_type}.merged.exp.kallisto.txt"
    output:
        "{cancer_type}.done"
    params:
        # cluster = "-l h_vmem=32G "
        cluster="--mem=32G -t 06:00:00"
    shell:
        """
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/gcc/9.2.0/lib64:/cm/local/apps/gcc/9.2.0/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R/lib:/cm/shared/apps_chop/R/4.0.2/lib64/R:/cm/shared/apps_chop/jags/4.3.0/lib:/cm/shared/apps_chop/gsl/2.5/lib:/cm/shared/apps_chop/gdal/2.4.0/lib:/cm/shared/apps_chop/hdf5/1.10.4/lib:/cm/shared/apps_chop/jdk/1.8.0_144/jre/lib/amd64/server:/cm/shared/apps_chop/jdk/1.8.0_144/lib:/cm/shared/apps_chop/perl/5.26.1/lib:/cm/shared/easybuild/software/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1/lib:/cm/shared/easybuild/software/FFTW/3.3.8-gompi-2018b/lib:/cm/shared/easybuild/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/lib:/cm/shared/easybuild/software/hwloc/1.11.10-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libpciaccess/0.14-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/libxml2/2.9.8-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/XZ/5.2.4-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/numactl/2.0.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/zlib/1.2.11-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/binutils/2.30-GCCcore-7.3.0/lib:/cm/shared/easybuild/software/GCCcore/7.3.0/lib/gcc/x86_64-pc-linux-gnu/7.3.0:/cm/shared/easybuild/software/GCCcore/7.3.0/lib64:/cm/shared/easybuild/software/GCCcore/7.3.0/lib:/cm/shared/apps_chop/sam-bcf-tools/1.6/lib:/cm/shared/apps/uge/8.6.12/lib/lx-amd64
        module load R
        module load GCC
        Rscript {REPO_DIR}/generate_SPARKS_object.R \
        --psi_SE rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.SE.csv \
        --psi_A3SS rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A3SS.csv \
        --psi_A5SS rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.A5SS.csv \
        --psi_RI rMATS_naive_matrix/filtered_psi.matrix.{wildcards.cancer_type}.RI.csv \
        --exp {input.exp_matrix} \
        --cancer_type {wildcards.cancer_type} \
        --exon_anno_SE rMATS_naive_matrix/AS_annotation.SE.{wildcards.cancer_type}.txt \
        --exon_anno_A3SS rMATS_naive_matrix/AS_annotation.A3SS.{wildcards.cancer_type}.txt \
        --exon_anno_A5SS rMATS_naive_matrix/AS_annotation.A5SS.{wildcards.cancer_type}.txt \
        --exon_anno_RI rMATS_naive_matrix/AS_annotation.RI.{wildcards.cancer_type}.txt

        touch {wildcards.cancer_type}.done

        """

