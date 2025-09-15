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

SAMPLE_INDEX = get_sampleids_from_path_pattern('Study_dataset1/*')


rule all:
	input:
        	expand("{home}/refined/{sampleID}",sampleID=SAMPLE_INDEX, home=HOME)

####-------------####
#### 5. Refine 
####-------------####


###############
## BinRefine ##
###############
# outputs are metawrap_50_10_bins.contig, metawrap_50_10_bins.stats, concoct.stats, maxbin.stats, metabat.stats,
# figures, metawrap_50_10_bins


rule binRefine:
    input:
        concoct = directory("{home}/concoct/{sampleID}/{sampleID}.concoct-bins"),
        metabat = directory("{home}/metabat/{sampleID}/{sampleID}.metabat-bins"),
        maxbin = directory("{home}/maxbin/{sampleID}/{sampleID}.maxbin-bins")
    output:
        directory('{home}/refined/{sampleID}')
    envmodules:
        "bioinfo-tools",
        "conda"
    benchmark:
        '{home}/benchmark/{sampleID}.binRefine.benchmark.txt'
    log:
        '{home}/logs/{sampleID}/{sampleID}_binRefine.log'
    shell:
        """
        echo -e "$(date)\nSection starts\n BinRefine \n"

        # Create output folder
        mkdir -p {output}
        idvar=$(echo $(basename {input.concoct})|sed 's/.concoct-bins//g')
        #mkdir -p {config[path][scratch]}/{config[folder][refined]}/{{idvar}}
        #cd {config[path][scratch]}/{config[folder][refined]}/{{idvar}}
        mkdir -p /cfs/klemming/scratch/x/xinmengl/refined/${{idvar}}
	cd /cfs/klemming/scratch/x/xinmengl/refined/${{idvar}}
	
        #set +e
        #set +u
        #conda activate --stack {config[envs][metawrap]}
        #set -e
        #set -u
        
        # Copy files to tmp
        echo "Copying bins from CONCOCT, metabat2, and maxbin2 ... "
        cp -r {input.concoct} {input.metabat} {input.maxbin} ./

        echo "Renaming bin folders to avoid errors with metaWRAP ... "
        rm -rf ${{idvar}}.concoct ${{idvar}}.metabat ${{idvar}}.maxbin
        mv $(basename {input.concoct}) $(echo $(basename {input.concoct})|sed 's/-bins//g') 
        mv $(basename {input.metabat}) $(echo $(basename {input.metabat})|sed 's/-bins//g') 
        mv $(basename {input.maxbin}) $(echo $(basename {input.maxbin})|sed 's/-bins//g')
        
        echo "Running metaWRAP bin refinement module ... "

        # take 1h for sample1 with 48 cores
        #set +e
        #set +u
        metaWRAP bin_refinement -o ./ \
            -A $(echo $(basename {input.concoct})|sed 's/-bins//g') \
            -B $(echo $(basename {input.metabat})|sed 's/-bins//g') \
            -C $(echo $(basename {input.maxbin})|sed 's/-bins//g') \
            -t {config[cores][refine]} \
            -m {config[params][refineMem]} \
            -c {config[params][refineComp]} \
            -x {config[params][refineCont]} \
            2> {log}
        #set -u
        #set -e
        
        rm -rf $(echo $(basename {input.concoct})|sed 's/-bins//g') $(echo $(basename {input.metabat})|sed 's/-bins//g') $(echo $(basename {input.maxbin})|sed 's/-bins//g') work_files
        mv * {output}

        #cd /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/metaGEM/workflow/scratch/refined_bins
        #metaWRAP bin_refinement -o ./ -A sample1.concoct -B sample1.maxbin -C sample1.metabat -t 48 -m 1600 -c 50 -x 10
        echo "Bin Refining Done at $(date)."
        """

