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

qfiltered_pattern = f"{root_path}/{config['folder']['qfiltered']}/*"
print(f"DEBUG: Looking for samples in: {qfiltered_pattern}")
SAMPLE_INDEX = get_sampleids_from_path_pattern(qfiltered_pattern)
print(f"DEBUG: Found samples: {SAMPLE_INDEX}")

rule all:
    input:
        expand(f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz", sampleID=SAMPLE_INDEX)



####-------------------------------------####
#### Megahit
####------------------------------------####
rule megahit:
    input:
        R1 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
    output:
        assembly = f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz"
    threads: config['cores']['megahit']
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.megahit.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_megahit.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Megahit \n"

        mkdir -p $(dirname {output.assembly})

        idvar=$(echo $(basename $(dirname {output.assembly})))
        ihome=$(echo $(dirname $(dirname {output.assembly})))
        echo -e "\nCreating temporary directory ${{ihome}}/${{idvar}} ..."
        mkdir -p ${{ihome}}/${{idvar}}
        cd ${{ihome}}/${{idvar}}

        echo -n "Copying qfiltered reads to ${{ihome}}/${{idvar}} ... "
        cp {input.R1} {input.R2} ./
        echo -n "Running megahit ..., assemble mode is {config[params][assemblyPreset]} ... "

        rm -rf tmp

        megahit -1 $(basename {input.R1}) \
            -2 $(basename {input.R2}) \
            -o tmp \
            -t {config[cores][megahit]} \
            --presets {config[params][assemblyPreset]} \
            --verbose \
            --min-contig-len {config[params][assemblyMin]}
        
        echo "done. "
        echo "Renaming assembly ... "
        mv tmp/final.contigs.fa contigs.fasta

        sed -i 's/ /-/g' contigs.fasta
        gzip contigs.fasta

        rm -rf tmp/

        echo "Megahit assembling Done at $(date). "
        """
