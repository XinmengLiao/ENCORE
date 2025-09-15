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

SAMPLE_INDEX = get_sampleids_from_path_pattern('Validation1/*')

rule all:
	input:
        	expand("{home}/concoct/{sampleID}/{sampleID}.concoct-bins", sampleID=SAMPLE_INDEX, home=HOME) ## Concoct


#############
## Concoct ##
#############
# scikit-learn v1.5 cannot be used with Concoct, need to change to v1.1.0
# outputs are 0-34.fa
# need to move all these files out of the folders: args.txt, clustering_gt1000.csv, clustering_merge.csv,
# original_data_gt1000.csv, PCA_components/transformed_data_gt1000.csv

rule concoct:
    input:
        contigs = "{home}/megahit_assembled/{sampleID}/contigs.fasta.gz"
    output:
        directory("{home}/concoct/{sampleID}/{sampleID}.concoct-bins")
    benchmark:
        "{home}/benchmark/{sampleID}.concoct.benchmark.txt"
    log:
        "{home}/logs/{sampleID}/{sampleID}_concoct.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Concoct \n"
        
        # Actiavte the conda environment 
        #set +e
        #set +u 
        #source activate metagem 
        #set -u 
        #set -e

        # Create output folder
        mkdir -p $(dirname {output})

        # Make job specific scratch dir
        sample=$(echo $(basename $(dirname {input.contigs})))
        echo -e "\nCreating temporary directory {config[path][scratch]}/{config[folder][concoct]}/${{sample}} ... "
        mkdir -p {config[path][scratch]}/{config[folder][concoct]}/${{sample}}

        # Move into scratch dir
        cd {config[path][scratch]}/{config[folder][concoct]}/${{sample}}

        # Copy files
        cp {input.contigs} $(dirname {output})/cov/coverage_table.tsv ./

        echo "Unzipping assembly ... "
        gunzip -f $(basename {input.contigs})

        echo "Adujusting coverage table."
        cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][adjust_coverage]} ./
        python3 adjust_coverage.py

        echo -e "Done. \nCutting up contigs (default 10kbp chunks) ... "
        cut_up_fasta.py -c {config[params][cutfasta]} -o 0 -m $(echo $(basename {input.contigs})|sed 's/.gz//') > assembly_c10k.fa
        
        echo -e "\nRunning CONCOCT ... "
        concoct --coverage_file $(dirname {output})/cov/coverage_table.tsv \
            --composition_file assembly_c10k.fa \
            -b $(basename $(dirname {output})) \
            -t {config[path][cores]}/{config[folder][concoct]} \
            -c 400 -i 350 \
            2> {log}
            
        echo -e "\nMerging clustering results into original contigs ... "
        merge_cutup_clustering.py $(basename $(dirname {output}))_clustering_gt1000.csv > $(basename $(dirname {output}))_clustering_merged.csv
        
        echo -e "\nExtracting bins ... "
        mkdir -p $(basename {output})
        extract_fasta_bins.py $(echo $(basename {input.contigs})|sed 's/.gz//') $(basename $(dirname {output}))_clustering_merged.csv --output_path $(basename {output})
		
        # Move final result files to output folder
        mkdir -p $(dirname {output})
        mv $(basename {output}) *.txt *.csv $(dirname {output})

	    echo "Concoct Done at $(date)."
        """
