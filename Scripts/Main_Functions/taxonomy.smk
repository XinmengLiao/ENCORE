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
		expand(f"{config['path']['root']}/{config['folder']['classification']}/{sampleID}", sampleID=SAMPLE_INDEX)


####------------------------------------####
#### Taxonomic assignment with GTDB-tk
####------------------------------------####


#########################
## GTDBtk for taxonomy ## (can run, but need to set the path and other parameters first)
#########################
rule GTDBtk:
    input:
        f"{config['path']['root']}/{config['folder']['dnaBins']}/{{sampleID}}"
    output:
        directory(f"{config['path']['root']}/{config['folder']['classification']}/{{sampleID}}")
    message:
        """
    # For ENCORE validation and application on PD cohorts, GTDB-Tk release220 was used. 
    # GTDB-Tk release220 full database contains ~100G of external data which needs to be downloaded manually and extracted.
    # Other versions of GTDB-Tk datasets can be found here: https://ecogenomics.github.io/GTDBTk/installing/index.html
    # PATH should be added to the environment variable: `GTDBTK_DATA_PATH="[PATH/TO/GTDBTK]/"`
        """
    shell:
        """
        echo -e "$(date)\nSection starts\n GTDBtk \n"  

        # Activate the env
        set +e;
        set +u;
        source metagem
        set -u;
        set -e

        # Make sure the output directory exists
        mkdir -p {output}

        #export GTDBTK_DATA_PATH={config[path][gtdbtk]}
	    set +u
        export PYTHONPATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/lib/python3.12/site-packages:$PYTHONPATH
	    export PATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/bin:$PATH
	    set -u
	    cd {input}

        gtdbtk classify_wf --genome_dir ./ \
            --out_dir GTDBtk \
            -x fa \
            --cpus {config[cores][gtdbtk]} \
            --mash_db ./ \
            2> {log}
        
        idvar=$(basename {output})
        for file in GTDBtk/*summary.tsv; do
            filename=$(basename "$file")
            cp $file {output}/${{idvar}}_${{filename}}
        done

        cd {output}
        # Combine all summary files
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' *.summary.tsv > all.tsv
        mv all.tsv ${{idvar}}_gtdbtk_summary.tsv

        echo "GTDBtk taxonomic assignment done at $(date). "
        """



