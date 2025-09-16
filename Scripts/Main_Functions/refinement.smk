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
        	expand(f"{config['path']['root']}/{config['folder']['refined']}/{sampleID}", sampleID=SAMPLE_INDEX)
                
####----------------####
#### Bin Refinement ####
####----------------####

rule binRefine:
    input:
        concoct = directory(f"{config['path']['root']}/{config['folder']['concoct']}/{{sampleID}}/{{sampleID}}.concoct-bins"),
        metabat = directory(f"{config['path']['root']}/{config['folder']['metabat']}/{{sampleID}}/{{sampleID}}.metabat-bins"),
        maxbin = directory(f"{config['path']['root']}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins")
    output:
        directory(f"{config['path']['root']}/{config['folder']['refined']}/{{sampleID}}")
    shell:
        """
        echo -e "$(date)\nSection starts\n BinRefine \n"

        # Create output folder
        mkdir -p {output}
        idvar=$(echo $(basename {input.concoct})|sed 's/.concoct-bins//g')
        mkdir -p {config[path][scratch]}/{config[folder][refined]}/${{idvar}}
        cd {config[path][scratch]}/{config[folder][refined]}/${{idvar}}
	
        set +e
        set +u
        conda activate --stack {config[envs][metawrap]}
        set -e
        set -u
        
        # Copy files to tmp
        echo "Copying bins from CONCOCT, metabat2, and maxbin2 ... "
        cp -r {input.concoct} {input.metabat} {input.maxbin} ./

        echo "Renaming bin folders to avoid errors with metaWRAP ... "
        rm -rf ${{idvar}}.concoct ${{idvar}}.metabat ${{idvar}}.maxbin
        mv $(basename {input.concoct}) $(echo $(basename {input.concoct})|sed 's/-bins//g') 
        mv $(basename {input.metabat}) $(echo $(basename {input.metabat})|sed 's/-bins//g') 
        mv $(basename {input.maxbin}) $(echo $(basename {input.maxbin})|sed 's/-bins//g')
        
        echo "Running metaWRAP bin refinement module ... "

        set +e
        set +u
        metaWRAP bin_refinement -o ./ \
            -A $(echo $(basename {input.concoct})|sed 's/-bins//g') \
            -B $(echo $(basename {input.metabat})|sed 's/-bins//g') \
            -C $(echo $(basename {input.maxbin})|sed 's/-bins//g') \
            -t {config[cores][refine]} \
            -m {config[params][refineMem]} \
            -c {config[params][refineComp]} \
            -x {config[params][refineCont]} \
            2> {log}
        set -u
        set -e
        
        rm -rf $(echo $(basename {input.concoct})|sed 's/-bins//g') $(echo $(basename {input.metabat})|sed 's/-bins//g') $(echo $(basename {input.maxbin})|sed 's/-bins//g') work_files
        mv * {output}

        echo "Bin Refining Done at $(date)."
        """

