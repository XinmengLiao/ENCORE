configfile: "ENCORE_config.yaml"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern(f"{config['path']['root']}/{config['folder']['data']}/*")

rule all:
    input:
        expand(f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX), 
        expand(f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX)



####--------####
#### Fastp
####--------####

rule fastp:
    input:
        R1 = f"{config['path']['root']}/{config['folder']['data']}/{{sampleID}}/{{sampleID}}_R1.fastq.gz",
        R2 = f"{config['path']['root']}/{config['folder']['data']}/{{sampleID}}/{{sampleID}}_R2.fastq.gz"
    output:
        R1 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
	shell:
		"""
        echo -e "$(date)\nSection starts\n ***** Fastp ***** \n"

        set +e
        set +u 
        source activate metagem 
        set -u 
        set -e

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

 ## Skip -- qfilterVis