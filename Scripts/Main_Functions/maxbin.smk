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
        expand(f"{config['path']['root']}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins", sampleID=SAMPLE_INDEX)



#################
## MaxbinCross ##
#################
rule maxbinCross:
    input:
        assembly = f"{config['path']['root']}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz",
        depth = f"{config['path']['root']}/{config['folder']['maxbin']}/{{sampleID}}/cov"
    output:
        directory(f"{config['path']['root']}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins")
    shell:
        """
        echo -e "$(date)\nSection starts\n Maxbin \n"
        
        # Create output folder
        mkdir -p $(dirname {output})

        # Make job specific scratch dir
        fsampleID=$(echo $(basename $(dirname {input.assembly})))
	
        mkdir -p {config['path']['scratch']}/{config['folder']['maxbin']}/${{fsampleID}}
        cd {config['path']['scratch']}/{config['folder']['maxbin']}/${{fsampleID}} 

        # Copy files to tmp
        cp -r {input.assembly} {input.depth}/*.depth ./

        echo -e "\nUnzipping assembly ... "
        gunzip -f $(basename {input.assembly})

        echo -e "\nGenerating list of depth files based on crossMapSeries rule output ... "
        find . -name "*.depth" > abund.list
        
        echo -e "\nRunning maxbin2 ... "
        run_MaxBin.pl -thread {config[cores][maxbin]} \
            -contig contigs.fasta \
            -out $(basename $(dirname {output})) \
            -abund_list abund.list \
            -prob_threshold 0.9 \
            -min_contig_length 1000 \
            -max_iteration 10
        
        rm *.depth abund.list contigs.fasta

        # Manage Data
        mkdir -p $(basename {output})
	    mv *.fasta $(basename {output})
        mv * $(dirname {output})

        echo "Maxbin Done at $(date). "
        """

