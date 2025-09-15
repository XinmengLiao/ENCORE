configfile: "GEMconfig.yaml"
#SAMPLE_INDEX = {"Feces_20004_Visit3_S9"}
HOME="/cfs/klemming/projects/supr/naiss2023-23-637/ad"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern('Validation1/*')

rule all:
	input:
        	expand("{home}/qfiltered/{sampleID}/{sampleID}_R1_qfiltered.fastq.gz", home = HOME, sampleID = SAMPLE_INDEX), ## Fastp
        	expand("{home}/qfiltered/{sampleID}/{sampleID}_R2_qfiltered.fastq.gz", home = HOME, sampleID = SAMPLE_INDEX)



####--------------------------------------####
#### 1. Fastp (no problem, run smoothly)
####--------------------------------------####

rule fastp:
	input:
		R1 = "{home}/dataset1/{sampleID}/{sampleID}_R1.fastq.gz",
		R2 = "{home}/dataset1/{sampleID}/{sampleID}_R2.fastq.gz"
	output:
		R1 = "{home}/qfiltered/{sampleID}/{sampleID}_R1_qfiltered.fastq.gz",
		R2 = "{home}/qfiltered/{sampleID}/{sampleID}_R2_qfiltered.fastq.gz"
	threads: 32
	benchmark:
		"{home}/benchmark/{sampleID}.fastp.benchmark.txt"
	shell:
		"""
        echo -e "$(date)\nSection starts\n ***** Fastp ***** \n"

        #module load bioinfo-tools
        #module load conda
        #metagem # activate conda metagem env

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