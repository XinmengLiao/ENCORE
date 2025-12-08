import os
import glob

configfile: os.path.join(workflow.basedir, "../../ENCORE_config.yaml")

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

DATA_FOLDER = config.get('DATA_FOLDER', 'Toy_Dataset')

# Build the root path - prefer path_root, then OUTPUT_DIR, then config value
root_path = config.get('path_root') or config.get('OUTPUT_DIR') or config['path']['root']

# Handle relative paths for root_path
if not os.path.isabs(root_path):
    root_path = os.path.join(os.path.dirname(workflow.basedir), root_path)

assemblies_pattern = f"{root_path}/{config['folder']['assemblies']}/*"
print(f"DEBUG: Looking for samples in: {assemblies_pattern}")
SAMPLE_INDEX = get_sampleids_from_path_pattern(assemblies_pattern)
print(f"DEBUG: Found samples: {SAMPLE_INDEX}")

rule all:
    input:
        expand(f"{root_path}/{config['folder']['concoct']}/{{sampleID}}/{{sampleID}}.concoct-bins", sampleID=SAMPLE_INDEX)


#############
## Concoct ##
#############
rule concoct:
    input:
        contigs = f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz",
        coverage = f"{root_path}/{config['folder']['concoct']}/{{sampleID}}/cov/coverage_table.tsv"
    output:
        directory(f"{root_path}/{config['folder']['concoct']}/{{sampleID}}/{{sampleID}}.concoct-bins")
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.concoct.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_concoct.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Concoct \n"

        sample=$(echo $(basename $(dirname {input.contigs})))
        echo -e "\nCreating temporary directory {root_path}/{config[folder][crossMap]}/${{sample}} ... "
        mkdir -p {root_path}/{config[folder][crossMap]}/${{sample}}

        cd {root_path}/{config[folder][crossMap]}/${{sample}}

        cp {input.contigs} {input.coverage} ./

        echo "Unzipping assembly ... "
        gunzip -f $(basename {input.contigs})

        echo "Adjusting coverage table."
        cp $(dirname {root_path})/{config[path][scripts]}/PythonScripts/adjust_coverage.py ./
        python3 adjust_coverage.py

        echo -e "Done. \nCutting up contigs (default 10kbp chunks) ... "
        cut_up_fasta.py -c {config[params][cutfasta]} -o 0 -m $(echo $(basename {input.contigs})|sed 's/.gz//') > assembly_c10k.fa

        echo -e "\nRunning CONCOCT ... "
        concoct --coverage_file coverage_table.tsv \
            --composition_file assembly_c10k.fa \
            -b ${{sample}} \
            -t {config[cores][concoct]} \
            -c {config[params][concoct]} \
            -i {config[params][concoctMinLength]}

        echo -e "\nMerging clustering results into original contigs ... "
        merge_cutup_clustering.py ${{sample}}_clustering_gt1000.csv > ${{sample}}_clustering_merged.csv

        echo -e "\nExtracting bins ... "
        mkdir -p ${{sample}}.concoct-bins
        extract_fasta_bins.py $(echo $(basename {input.contigs})|sed 's/.gz//') ${{sample}}_clustering_merged.csv --output_path ${{sample}}.concoct-bins

        mkdir -p {output}
        mv ${{sample}}.concoct-bins/* {output}/

        echo "Concoct Done at $(date)."
        """