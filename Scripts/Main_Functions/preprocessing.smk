"""
ENCORE Preprocessing Pipeline
Integrated workflow: FastQC -> fastp -> MEGAHIT Assembly
"""

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

# Handle relative paths for root_path
if not os.path.isabs(root_path):
    root_path = os.path.join(os.path.dirname(workflow.basedir), root_path)

raw_pattern = f"{root_path}/{DATA_FOLDER}/*"
print(f"DEBUG: Looking for samples in: {raw_pattern}")
SAMPLE_INDEX = get_sampleids_from_path_pattern(raw_pattern)
print(f"DEBUG: Found samples: {SAMPLE_INDEX}")

# ================================
# Include individual workflow modules
# ================================
include: "fastqc.smk"
include: "fastp.smk"
include: "megahit.smk"

# ================================
# Pipeline-level rule
# ================================
rule all:
    input:
        # FastQC outputs
        expand(f"{root_path}/{config['folder']['quality']}/{{sampleID}}/{{sampleID}}_R1_fastqc.zip", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['quality']}/{{sampleID}}/{{sampleID}}_R2_fastqc.zip", sampleID=SAMPLE_INDEX),
        # Fastp outputs
        expand(f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz", sampleID=SAMPLE_INDEX),
        # MEGAHIT assembly outputs
        expand(f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz", sampleID=SAMPLE_INDEX)
