configfile: "GEMconfig.yaml"
HOME="/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/metaGEM/workflow"


import os
import glob



def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern('Study_dataset/*')


rule all:
	input:
		expand("{home}/GEMs/{sampleID}", home=HOME, sampleID=SAMPLE_INDEX) # CarveMe



####------------------------####
#### 6. GEMs Constructions
####------------------------####


#############
## Carveme ## (can run, but need to set other parameters first)
#############
# outputs are sample.bin.1-4.s.xml, sample.bin.1-4.p.xml
# need to pip install carveme and conda install diamond (sequence aligner)
# Here all idvar should not use {}(try if works). But here did not remove, can try again. 
# every time when start to use Carveme, need to export the path:
# export PYTHONPATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/lib/python3.12/site-packages:$PYTHONPATH
# export PATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/bin:$PATH
# if can not find carveme, then `pip show carveme`
# if does ncarveot work, then just uninstall carveme and install again. 

rule carveme:
    input:
        bin = "{home}/protein_bins/{sampleID}",
        media = "{home}/scripts/media_db.tsv"
    output:
        directory('{home}/GEMs/{sampleID}')
    benchmark:
        "{home}/benchmark/{sampleID}/{sampleID}.carveme.benchmark.txt"
    log:
        "{home}/logs/{sampleID}/{sampleID}_carveme.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Carveme \n"
        
        #set +e
        #set +u
        #source activate metagem
        #set -u
        #set -e

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

