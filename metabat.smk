configfile: "GEMconfig.yaml"
HOME = "/cfs/klemming/projects/supr/naiss2023-23-637/ecgem/metaGEM/workflow"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern(f'{HOME}/dataset[1-4]/*')

rule all:
	input:
            expand("{home}/metabat/{sampleID}/{sampleID}.metabat-bins", home=HOME, sampleID=SAMPLE_INDEX), # Metabat output


##################
## MetabatCross ##
##################
# output is sample.1-10.fa

rule metabatCross:
    input:
        assembly = "{home}/megahit_assembled/{sampleID}/contigs.fasta.gz",
        depth = "{home}/metabat/{sampleID}/cov"
    output:
        directory("{home}/metabat/{sampleID}/{sampleID}.metabat-bins")
    benchmark:
        "{home}/benchmark/{sampleID}.metabat.benchmark.txt"
    log:
        "{home}/logs/{sampleID}/{sampleID}_metabatcross.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Metabat \n"

        # Actiavte the conda environment 
        #set +e
        #set +u 
        #source activate metagem 
        #set -u 
        #set -e
        
        # Make job specific scratch dir
        echo -e "\nCreating temporary directory {config[path][scratch]}/{config[folder][metabat]}/{wildcards.sampleID} ... "
        mkdir -p {config[path][scratch]}/{config[folder][metabat]}/{wildcards.sampleID}

        # Move into scratch dir
        echo -e "\nMoving into temporary directory {config[path][scratch]}/{config[folder][metabat]}/{wildcards.sampleID} ... "
        cd {config[path][scratch]}/{config[folder][metabat]}/{wildcards.sampleID}

        # Create output folder
        mkdir -p {output}
        
        # Copy with verification
        cp {input.assembly} ./contigs.fasta.gz
        ls -l ./contigs.fasta.gz >> {log}
        
        cp {input.depth}/*all.depth ./
        ls -l ./*.all.depth >> {log}
        
        # Unzip with verification
        gunzip -f ./contigs.fasta.gz
        ls -l ./contigs.fasta >> {log}
        
        # Check file sizes
        du -h ./contigs.fasta >> {log}
        du -h ./*.all.depth >> {log}
        
        # Run metabat2 with memory info
        /usr/bin/time -v metabat2 -i contigs.fasta \
            -a {wildcards.sampleID}.all.depth \
            -s {config[params][metabatMin]} \
            -v \
            --seed {config[params][seed]} \
            -t 0 \
            -m {config[params][minBin]} \
            -o {wildcards.sampleID} \
            2>> {log}

        # Check and move .fa files if they exist
        if compgen -G "*.fa" > /dev/null; then
            mv *.fa {output}
            echo "Moved .fa files to output directory" >> {log}
        else
            mkdir -p {output}
            echo "No .fa files found to move" >> {log}
        fi

        echo "Metabat done at $(date). "
        """
