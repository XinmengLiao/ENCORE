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
        expand(f"{config['path']['root']}/{config['folder']['GEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX) # CarveMe



####------------------------####
#### 6. GEMs Constructions
####------------------------####


#############
## Carveme ## (can run, but need to set other parameters first)
#############
rule carveme:
    input:
        bin = f"{config['path']['root']}/{config['folder']['proteinBins']}/{{sampleID}}",
        media = f"{config['path']['root']}/{config['folder']['scripts']}/{config['scripts']['carveme']}"
    output:
        directory(f"{config['path']['root']}/{config['folder']['GEMs']}/{{sampleID}}")
    shell:
        """
        echo -e "$(date)\nSection starts\n Carveme \n"
        
        set +e
        set +u
        source activate metagem
        set -u
        set -e

        # Make sure output folder exists
        idvar=$(echo $(basename {output}))
        mkdir -p {output}
        mkdir -p {config[path][scratch]}/{config[folder][GEMs]}/${{idvar}}
        cd {config[path][scratch]}/{config[folder][GEMs]}/${{idvar}}

        # Copy files
        cp {input.media} {input.bin}/* ./
        
        # Add Carveme to PATH

        for name in *.faa; do
            binID=$(echo $name | sed 's/.faa//g')
            echo "Begin carving GEM for sample: ${{idvar}} protein bin: ${{binID}} ... "
            carve $name \
                -g {config[params][carveMedia]} \
                -v \
                --mediadb $(basename {input.media}) \
                --fbc2 \
                -o ${{binID}}.xml \
                2> {log}
            echo "${{idvar}} protein bin: ${{binID}} finished constructing GEMs at $(date). ";
        done

	    mkdir -p {output}
        mv *.xml {output}
        echo "Done carving GEM (Carveme) for all protein bins of Sample: ${{idvar}}. Now start to generate the statistic for GEMs.  "
        

        ## GEMs statistics
        cd {output}
        while read model;do 
            id=$(echo $(basename $model)|sed 's/.xml//g'); 
            mets=$(less $model| grep "species metaid="|cut -d ' ' -f 8|sort|uniq|wc -l);
            rxns=$(less $model|grep -c 'reaction metaid=');
            genes=$(less $model|grep -E 'fbc:geneProduct.*fbc:id='|wc -l);
            echo "Model: $id has $mets mets, $rxns reactions, and $genes genes ... "
            echo "$id $mets $rxns $genes" >> GEMs.stats;
        done< <(find . -name "*.xml")
	
	    mkdir -p {config[path][root]}/{config[folder][stats]}/${{idvar}}
        mv GEMs.stats {config[path][root]}/{config[folder][stats]}/${{idvar}}
        cd {config[path][root]}/{config[folder][stats]}/${{idvar}}
        Rscript {config[path][root]}/{config[folder][scripts]}/{config[scripts][modelVis]}
        #rm Rplots.pdf # Delete redundant pdf file
        echo "GEMs statistic checking done at $(date). "

        """

