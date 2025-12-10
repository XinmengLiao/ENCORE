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
        f"{root_path}/{config['folder']['network']}/MMNetwork.txt"


####-------------------------------------------####
#### Generation of Microbe-Metabolite Network
####-------------------------------------------####

rule MMNetwork:
    input:
        dataDir =  f"{root_path}/{config['folder']['assemblies']}",
        ecgemDir =  f"{root_path}/{config['folder']['ecGEMs']}", 
        taxoDir = f"{root_path}/{config['folder']['classification']}"
    output:
        f"{root_path}/{config['folder']['network']}/MMNetwork.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        cd $(dirname {output})

        # Process each sample
        for sampledir in {input.dataDir}/*; do
            sampleid=$(basename "$sampledir")
            echo "Processing sample: $sampleid"

            mkdir -p $(dirname {output})/taxonomy
            mkdir -p $(dirname {output})/fva

            # Copy all taxonomy files for this sample
            cp {input.taxoDir}/${{sampleid}}/*_gtdbtk_summary.tsv ./taxonomy/

            # Process each bin in this sample
            for bindir in {input.ecgemDir}/${{sampleid}}/bin* ; do
                if [ -d "$bindir" ]; then                    
                    cp $bindir/output/*FVA.txt ./fva
                fi
            done
        done

        echo "Now creating Microbe-Metabolite Network..."
        Rscript $(dirname {root_path})/{config[Rscripts][createMMNetwork]} $(dirname {output})/taxonomy $(dirname {output})/fva $(dirname {output})

        """