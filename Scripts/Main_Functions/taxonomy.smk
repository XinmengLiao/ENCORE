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
        expand(f"{root_path}/{config['folder']['classification']}/{{sampleID}}", sampleID=SAMPLE_INDEX)


####------------------------------------####
#### Taxonomic assignment with GTDB-tk
####------------------------------------####


#########################
## GTDBtk for taxonomy ## 
#########################
rule GTDBtk:
    input:
        f"{root_path}/{config['folder']['dnaBins']}/{{sampleID}}"
    output:
        directory(f"{root_path}/{config['folder']['classification']}/{{sampleID}}")
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.taxonomy.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_taxonomy.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n GTDBtk \n"

        mkdir -p {output}

        # need to set the path and other parameters before running gtdbtk
        #export PYTHONPATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/lib/python3.12/site-packages:$PYTHONPATH
        #export PATH=/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/bin:$PATH

        cd {input}

        gtdbtk classify_wf --genome_dir ./ \
            --out_dir GTDBtk \
            -x fa \
            --cpus {config[cores][gtdbtk]} \
            --mash_db ./

        idvar=$(basename {output})
        for file in GTDBtk/*summary.tsv; do
            filename=$(basename "$file")
            cp $file {output}/${{idvar}}_${{filename}}
        done

        cd {output}
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' *.summary.tsv > all.tsv
        mv all.tsv ${{idvar}}_gtdbtk_summary.tsv

        echo "GTDBtk taxonomic assignment done at $(date)."
        """



