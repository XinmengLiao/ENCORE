import os
import glob

configfile: os.path.join(workflow.basedir, "../../ENCORE_config.yaml")

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

DATA_FOLDER = config.get('DATA_FOLDER', 'Toy_Dataset')

# Build the root path - prefer path_root, then OUTPUT_DIR, then config value
root_path = config.get('path_root') or config.get('OUTPUT_DIR') or config['path']['root']

# Handle relative paths
if not os.path.isabs(root_path):
    root_path = os.path.join(os.path.dirname(workflow.basedir), root_path)

SAMPLE_INDEX = get_sampleids_from_path_pattern(f"{root_path}/{DATA_FOLDER}/*")

rule all:
    input:
        expand(f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX), 
        expand(f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX)



####--------####
#### Fastp
####--------####

rule fastp:
    input:
        R1 = f"{root_path}/{DATA_FOLDER}/{{sampleID}}/{{sampleID}}_R1.fastq.gz",
        R2 = f"{root_path}/{DATA_FOLDER}/{{sampleID}}/{{sampleID}}_R2.fastq.gz"
    output:
        R1 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
    shell:
        """
        echo -e "$(date)\nSection starts\n ***** Fastp ***** \n"

        mkdir -p $(dirname {output.R1})
        fastp -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
            --thread {config[cores][fastp]} \
            -j $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).json \
            -h $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).html

        echo "Fastp quality filtering done at $(date). "
        """
