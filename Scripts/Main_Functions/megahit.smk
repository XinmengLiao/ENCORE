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
        expand(f"{config['path']['root']}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz", sampleID=SAMPLE_INDEX)



####-------------------------------------####
#### Megahit (no problem, run smoothly)
####------------------------------------####
rule megahit:
    input:
        R1 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
    output:
        f"{config['path']['root']}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz"
    threads: {config[cores][megahit]}
    benchmark:
        f"{config['path']['root']}/{config['folder']['benchmarks']}/{{sampleID}}.megahit.benchmark.txt"
    log:
        f"{config['path']['root']}/{config['folder']['logs']}/{{sampleID}}_megahit.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Megahit \n"
        
        # Actiavte the conda environment 
        set +e
        set +u 
        source activate metagem 
        set -u 
        set -e

        mkdir -p $(dirname {output})

        #Create the temporary folder
        idvar=$(echo $(basename $(dirname {output})))
        ihome=$(echo $(dirname $(dirname {output})))
        echo -e "\nCreating temporary directory ${{ihome}}/${{idvar}} ..."
        mkdir -p ${{ihome}}/${{idvar}}
        cd ${{ihome}}/${{idvar}}

        echo -n "Copying qfiltered reads to ${{ihome}}/${{idvar}} ... "
        cp {input.R1} {input.R2} ./
        echo -n "Running megahit ..., assemble mode is {config[params][assemblyPreset]} ... "

        # make sure there is not previous intermediate file 
        rm -rf tmp 

        megahit -1 $(basename {input.R1}) \
            -2 $(basename {input.R2}) \
            -o tmp \
            -t {config[cores][megahit]} \
            --presets {config[params][assemblyPreset]} \
            --verbose \
            --min-contig-len {config[params][assemblyMin]} \
            #--k-min  \
            #--k-max 99 \
            #--k-step 20 \
            2> {log}
        echo "done. "
        echo "Renaming assembly ... "
        mv tmp/final.contigs.fa contigs.fasta

        # Remove spaces from the contig headers and replace with hyphens
        sed -i 's/ /-/g' contigs.fasta
        gzip contigs.fasta
	
        rm -rf tmp/
        
        echo "Megahit assembling Done at $(date). "
	"""
