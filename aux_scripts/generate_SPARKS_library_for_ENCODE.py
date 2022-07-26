import pandas as pd
import os
import sys

snakedir = "/home/yangt3/xinglab/tools/SPARKS/"

input_df = pd.read_csv('metadata.tsv', sep='\t')

# define experiment type
experiment_type = "CRISPR"
# experiment_type = "shRNA"

# trim data
input_df = input_df.assign(Exp_type=experiment_type)
targets = [x.split('-')[0] for x in input_df['Experiment target']]
input_df = input_df.assign(Target=targets)

# aggregate replicate data
replicates = ["rep%s_%s" % (x['Biological replicate(s)'], x['Paired end']) for _, x in input_df.iterrows()]
input_df = input_df.assign(Replicate=replicates)

# determine control samples
control_list = []
for _, entry in input_df.iterrows():
    controls = entry['Controlled by'].split(',')
    if len(controls) > 1:
        replicate = entry['Biological replicate(s)']
        if str(replicate) == "1":
            control = controls[0]
        elif str(replicate) == "2":
            control = controls[1]
        else:
            sys.exit(333)
    else:
        control = controls[0]
    # trim control
    control_trimmed = control.split('/')[2]

    # add it to the list for attachment later
    control_list.append(control_trimmed)
input_df = input_df.assign(Control=control_list)

# merge all information for processing
sample_name = ["%s_%s_%s.%s" % (x['Biosample term name'],
                                x['Target'],
                                x['Exp_type'],
                                x['Replicate']) for _, x in input_df.iterrows()]

sample_project = ["%s_%s_%s" % (x['Biosample term name'],
                                x['Target'],
                                x['Exp_type']) for _, x in input_df.iterrows()]
input_df = input_df.assign(sample_file=sample_name)
input_df = input_df.assign(sample_project=sample_project)

control_dir = "/home/yangt3/xinglab/ENCODE/%s/control/" % experiment_type

# make script dir 
script_dir = "./scripts"
if not os.path.exists(script_dir):
    os.mkdir(script_dir)


# generate scripts for each file
projects = set(input_df['sample_project'].tolist())
for project in projects:
    command = ""

    command += "mkdir -p ../%s" % project + " \n "
    command += "mkdir -p ../%s/fastq" % project + " \n "
    command += "cd ../%s" % project + " \n "
    project_sample_df = input_df[input_df['sample_project'] == project]
    for _, entry in project_sample_df.iterrows():
        # download file
        sample_file = entry['sample_file']
        download_address = entry['File download URL']
        command += "curl -j -L -o ../%s/fastq/%s.fastq.gz %s \n " % (project, sample_file, download_address)

        # link control files
        control_file = entry['Control']
        command += "cp -v %s/%s.fastq.gz fastq/%s_Control_%s.%s.fastq.gz \n " % (control_dir, control_file, entry['Biosample term name'], entry['Exp_type'], entry['Replicate'])

    # copy snakemake framework file
    command += "ln -s %s/master.smk ./snakefile \n " % snakedir
    command += "ln -s %s/generate_SPARKS_object.R . \n " % snakedir
    
    # write configuration file
    command += "echo \"paired_end : True\" > ./config.yaml \n "
    command += "echo \"study : \\\"%s\\\" \" >> ./config.yaml \n " % project
    command += "echo \"read_length : \\\"100\\\" \" >> ./config.yaml \n " 
    command += "echo \"stat_mode : True\" >> ./config.yaml \n "
    command += "echo \"clip_mode : True\" >> ./config.yaml \n "
    command += "echo \"sample_definition:\" >> ./config.yaml \n "
    command += "echo \"    case : [\\\"%s\\\"]\" >> ./config.yaml \n " % project
    command += "echo \"    control : [\\\"%s_Control_%s\\\"]\" >> ./config.yaml \n " % (project.split("_")[0], project.split("_")[2])




    # generate snakemake commands
    command += "echo \"snakemake --cluster 'sbatch {params.cluster} -e ./logs/{rulename}.{jobid}.e%%A -o ./logs/{" \
               "rulename}.{jobid}.o%%A ' %s.done --jobs 30 --latency-wait 60 --jobname '%s.{rulename}.{jobid}' " \
               "--ri --keep-going\" > ../%s/run_START.%s.sh" % (project, project, project, project) + " \n "

    # generate new screen for snakemake
    command += "cd ../%s" % project + " \n "
    command += "screen -dm -S %s sh -c './run_START.%s.sh; exec bash'" % (project, project)

    # write command
    script_file = "%s/prepare_START_%s.sh" % (script_dir, project)
    with open(script_file, 'w') as f:
        f.write(command)

    print("Generated Scripts for %s" % project)
