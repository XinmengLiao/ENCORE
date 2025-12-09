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
        expand(f"{root_path}/{config['folder']['reassembled']}/{{sampleID}}", sampleID=SAMPLE_INDEX)

####################
## BinReassemble 
####################
rule binReassemble:
    input:
        R1 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz",
        refinedBins = directory(f"{root_path}/{config['folder']['refined']}/{{sampleID}}")
    output:
        reassembled_bin = directory(f"{root_path}/{config['folder']['reassembled']}/{{sampleID}}")
    params:
        scratch_dir = f"{root_path}/{config['path']['scratch']}/{config['folder']['reassembled']}"
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.reassembly.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_reassembly.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n BinReassemble \n"

        idvar=$(basename {output.reassembled_bin})
        mkdir -p {output.reassembled_bin}
        mkdir -p {params.scratch_dir}/${{idvar}}
        cd {params.scratch_dir}/${{idvar}}

        cp -r {input.refinedBins}/metawrap_*_bins {input.R1} {input.R2} ./

        rm -rf reassemblies/
        rm -rf reassembled_bins/
        rm -rf ${{idvar}}/

        echo "Running metaWRAP bin reassembly ... "
        metaWRAP reassemble_bins -o ${{idvar}} \
            -b metawrap_*_bins \
            -1 $(basename {input.R1}) \
            -2 $(basename {input.R2}) \
            -t {config[cores][reassemble]} \
            -m {config[params][reassembleMem]} \
            -c {config[params][reassembleComp]} \
            -x {config[params][reassembleCont]} \
            2>> {log}

        rm -rf metawrap_*_bins
        rm *.fastq.gz
        mv * {output.reassembled_bin}

        echo "BinReassembling done at $(date). Now start to extract protein bins."
        """