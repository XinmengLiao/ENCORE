configfile: "GEMconfig.yaml"
HOME="/cfs/klemming/projects/supr/naiss2023-23-637/ecgem/metaGEM/workflow"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern('dataset2/*')

rule all:
	input:
        	expand("{home}/maxbin/{sampleID}/{sampleID}.maxbin-bins", sampleID=SAMPLE_INDEX, home=HOME) ## Maxbin



#################
## MaxbinCross ##
#################
# output is sample.maxbin-bins.001-005.fasta

rule maxbinCross:
    input:
        assembly = "{home}/megahit_assembled/{sampleID}/contigs.fasta.gz",
        depth = "{home}/maxbin/{sampleID}/cov"
    output:
        directory("{home}/maxbin/{sampleID}/{sampleID}.maxbin-bins")
    benchmark:
        "{home}/benchmark/{sampleID}.maxbin.benchmark.txt"
    log:
        "{home}/logs/{sampleID}/{sampleID}_maxbincross.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Maxbin \n"
        
        # Create output folder
        mkdir -p $(dirname {output})

        # Make job specific scratch dir
        fsampleID=$(echo $(basename $(dirname {input.assembly})))
        #mkdir -p {config[path][scratch]}/{config[folder][maxbin]}/${{fsampleID}}
        #cd {config[path][scratch]}/{config[folder][maxbin]}/${{fsampleID}}
	
	mkdir -p /cfs/klemming/scratch/x/xinmengl/maxbin/${{fsampleID}}
	cd /cfs/klemming/scratch/x/xinmengl/maxbin/${{fsampleID}}	

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

