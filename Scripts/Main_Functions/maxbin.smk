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
        expand(f"{root_path}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins", sampleID=SAMPLE_INDEX)



#################
## MaxbinCross ##
#################
rule maxbinCross:
    input:
        assembly = f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz",
        depth = f"{root_path}/{config['folder']['maxbin']}/{{sampleID}}/cov"
    output:
        directory(f"{root_path}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins")
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.maxbin.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_maxbin.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Maxbin \n"

        mkdir -p $(dirname {output})

        fsampleID=$(echo $(basename $(dirname {input.assembly})))

        mkdir -p {root_path}/{config[path][scratch]}/maxbin/${{fsampleID}}
        cd {root_path}/{config[path][scratch]}/maxbin/${{fsampleID}}

        cp -r {input.assembly} {input.depth}/*.depth ./

        echo -e "\nUnzipping assembly ... "
        gunzip -f $(basename {input.assembly})

        echo -e "\nGenerating list of depth files based on crossMapSeries rule output ... "
        find . -name "*.depth" > abund.list

        echo -e "\nRunning maxbin2 ... "
        run_MaxBin.pl -thread {config[cores][maxbin]} \
            -contig contigs.fasta \
            -out $(basename $(dirname {output})) \
            -abund_list abund.list \
            -prob_threshold {config[params][maxbinProb]} \
            -min_contig_length {config[params][maxbinMinLength]} \
            -max_iteration {config[params][maxbinIteration]}

        rm *.depth abund.list contigs.fasta

        mkdir -p $(basename {output})
        mv *.fasta $(basename {output})
        mv * $(dirname {output})

        rm -rf {root_path}/{config[path][scratch]}/maxbin/${{fsampleID}}

        echo "Maxbin Done at $(date)."
        """

