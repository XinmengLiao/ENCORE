configfile: "GEMconfig.yaml"
SAMPLE_INDEX = {"Feces_20074_Visit3_S25"}
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
		expand("{home}/GTDBtk/{sampleID}", home=HOME, sampleID=SAMPLE_INDEX)


####---------------------------------------####
#### 7. Taxonomic assignment with GTDB-tk
####---------------------------------------####


#########################
## GTDBtk for taxonomy ## (can run, but need to set the path and other parameters first)
#########################
# GTDBtk can directly run in the metagem env
# but need to manually set the database path first: 
# `GTDBTK_DATA_PATH="/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/release220/"`
# numpy need the v1.23.1 version. `conda install -c conda-forge numpy=1.23.1`
# different version GTDBtk and the resource is imcompatible. release220 need a >=2.4 version
# installed gtdbtk v2.4 via pip in PDC. So can use 100G resource, but need to reinstall 'skani'
# In case of mess up, need to ensure other version of GTDBtk is not installed in conda. 


    # Previous METAGEM author information 
    # GTDB-Tk v2.1.1 requires ~63G of external data which needs to be downloaded and extracted. This can be done automatically, or manually.

    # Automatic:

    #     1. Run the command "download-db.sh" to automatically download and extract to:
    #         /home/xmliao/.conda/envs/gtdbtk/share/gtdbtk-2.1.1/db/

    # Manual:

    #     1. Manually download the latest reference data:
    #         wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz

    #     2. Extract the archive to a target directory:
    #         tar -xvzf gtdbtk_r207_v2_data.tar.gz -C "/path/to/target/db" --strip 1 > /dev/null
    #         rm gtdbtk_r207_v2_data.tar.gz

    #     3. Set the GTDBTK_DATA_PATH environment variable by running:
    #         conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"

rule GTDBtk:
    input:
        "{home}/dna_bins/{sampleID}"
    output:
        directory('{home}/GTDBtk/{sampleID}')
    benchmark:
        '{home}/benchmark/{sampleID}.GTDBtk.benchmark.txt'
    log:
        '{home}/logs/{sampleID}.GTDBtk.log'
    message:
        """

    Update 2024-09-05
    GTDB-TK v2.4 requires ~100G of external data which needs to be downloaded and extracted. This can be done automatically, or manually. 
    For more information, please refer to https://ecogenomics.github.io/GTDBTk/installing/index.html

    In PDC, v2.4 is installed by `pip install gtdbtk`, as the latest version currently is v2.4
    For other versions and the supported datasets, please check: https://ecogenomics.github.io/GTDBTk/installing/index.html

    Every time before using, if can not find gtdbtk in the system paths, then need to manually export the path by:
    `export PYTHONPATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/lib/python3.12/site-packages:$PYTHONPATH`
    `export PATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/bin:$PATH`

    Before running the GTDB-TK, need to set the GTDBTK_DATA_PATH environment variable: 
    `GTDBTK_DATA_PATH="/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/release220/"`

    
    Some dependencies need to be installed for running:
    (if using pip install gtdbtk, then all the dependencies need to manually install. Recommend using mamba install)
    1. numpy v1.23.1
    2. skani
    3. prodigal
    4. mash 

        """
    shell:
        """
        echo -e "$(date)\nSection starts\n GTDBtk \n"  

        # Activate the env
        #set +e;
        #set +u;
        #source activate gtdbtk-tmp;
        #set -u;
        #set -e

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



