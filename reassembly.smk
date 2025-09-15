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

SAMPLE_INDEX = get_sampleids_from_path_pattern('Study/*')


rule all:
	input:
        	expand("{home}/reassembled_bins/{sampleID}",sampleID=SAMPLE_INDEX, home=HOME) ## Bin Reassemble


####################
## BinReassemble 
####################
# need to use metawrap
# input refinedBins is metawrap_50_10_bins
# outputs are original_bins, original_bins.stats, reassembled_bins, reassembled_bins.checkm, 
# reassembled_bins.png, reassembled_bins.stats, reassembly_results.eps, reassembly_results.png

rule binReassemble:
    input:
        R1 = "{home}/qfiltered/{sampleID}/{sampleID}_R1_qfiltered.fastq.gz",
	R2 = "{home}/qfiltered/{sampleID}/{sampleID}_R2_qfiltered.fastq.gz",
        refinedBins = directory("{home}/refined/{sampleID}")
    output:
        reassembled_bin = directory('{home}/reassembled_bins/{sampleID}'),
    envmodules:
        "bioinfo-tools",
        "conda"
    benchmark:
        '{home}/benchmark/{sampleID}.binReassemble.benchmark.txt'
    log:
        "{home}/logs/{sampleID}/{sampleID}_binreassemble.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n BinReassemble \n"

        #set +e
        #set +u
        #source activate {config[envs][metagem]}
        #source activate {config[envs][metawrap]}
        #set -u
        #set -e

        idvar=$(basename {output.reassembled_bin})
        mkdir -p {output.reassembled_bin}
        #mkdir -p {config[path][scratch]}/{config[folder][reassembled]}/${{idvar}}
        mkdir -p /cfs/klemming/scratch/x/xinmengl/refined/${{idvar}}
        cd /cfs/klemming/scratch/x/xinmengl/refined/${{idvar}}

        #cd {config[path][scratch]}/{config[folder][reassembled]}/${{idvar}}
        cp -r {input.refinedBins}/metawrap_*_bins {input.R1} {input.R2} ./
        
        rm -rf reassemblies/
        rm -rf reassembled_bins/
        rm -rf ${{idvar}}/

        echo "Running metaWRAP bin reassembly ... "
        #set +e
        #set +u
        #conda activate --stack {config[envs][metawrap]}
        
        metaWRAP reassemble_bins -o ${{idvar}} \
            -b metawrap_*_bins \
            -1 $(basename {input.R1}) \
            -2 $(basename {input.R2}) \
            -t {config[cores][reassemble]} \
            -m {config[params][reassembleMem]} \
            -c {config[params][reassembleComp]} \
            -x {config[params][reassembleCont]} \
            2> {log}
        
        #set -u
        #set -e

        rm -rf metawrap_*_bins # remove the input 
        rm *.fastq.gz # remove the input 
        #rm -rf $(basename {output})/work_files
        mv * {output.reassembled_bin}
        
        echo "BinReassembling done at $(date). Now start to extrat protein bins. "

        """