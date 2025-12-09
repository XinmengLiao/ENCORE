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
        expand(f"{root_path}/{config['folder']['refined']}/{{sampleID}}", sampleID=SAMPLE_INDEX)
                
####----------------####
#### Bin Refinement ####
####----------------####

rule binRefine:
    input:
        concoct = directory(f"{root_path}/{config['folder']['concoct']}/{{sampleID}}/{{sampleID}}.concoct-bins"),
        metabat = directory(f"{root_path}/{config['folder']['metabat']}/{{sampleID}}/{{sampleID}}.metabat-bins"),
        maxbin = directory(f"{root_path}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins")
    output:
        directory(f"{root_path}/{config['folder']['refined']}/{{sampleID}}")
    params:
        scratch_dir = f"{root_path}/{config['path']['scratch']}/{config['folder']['refined']}"
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.refinement.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_refinement.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n BinRefine \n"

        mkdir -p {output}
        idvar=$(echo $(basename {input.concoct})|sed 's/.concoct-bins//g')
        mkdir -p {params.scratch_dir}/${{idvar}}
        cd {params.scratch_dir}/${{idvar}}

        echo "Copying bins from CONCOCT, metabat2, and maxbin2 ... "
        cp -r {input.concoct} {input.metabat} {input.maxbin} ./

        echo "Renaming bin folders to avoid errors with metaWRAP ... "
        rm -rf ${{idvar}}.concoct ${{idvar}}.metabat ${{idvar}}.maxbin
        mv $(basename {input.concoct}) $(echo $(basename {input.concoct})|sed 's/-bins//g')
        mv $(basename {input.metabat}) $(echo $(basename {input.metabat})|sed 's/-bins//g')
        mv $(basename {input.maxbin}) $(echo $(basename {input.maxbin})|sed 's/-bins//g')

        echo "Running metaWRAP bin refinement module ... "
        metaWRAP bin_refinement -o ./ \
            -A $(echo $(basename {input.concoct})|sed 's/-bins//g') \
            -B $(echo $(basename {input.metabat})|sed 's/-bins//g') \
            -C $(echo $(basename {input.maxbin})|sed 's/-bins//g') \
            -t {config[cores][refine]} \
            -m {config[params][refineMem]} \
            -c {config[params][refineComp]} \
            -x {config[params][refineCont]}

        rm -rf $(echo $(basename {input.concoct})|sed 's/-bins//g') $(echo $(basename {input.metabat})|sed 's/-bins//g') $(echo $(basename {input.maxbin})|sed 's/-bins//g') work_files
        mv * {output}

        echo "Bin Refining Done at $(date)."
        """

