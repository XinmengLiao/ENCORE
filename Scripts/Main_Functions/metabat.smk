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
        expand(f"{root_path}/{config['folder']['metabat']}/{{sampleID}}/{{sampleID}}.metabat-bins", sampleID=SAMPLE_INDEX)


##################
## MetabatCross ##
##################
rule metabatCross:
    input:
        assembly = f"{root_path}/{config['folder']['assemblies']}/{{sampleID}}/contigs.fasta.gz",
        depth = f"{root_path}/{config['folder']['metabat']}/{{sampleID}}/cov"
    output:
        directory(f"{root_path}/{config['folder']['metabat']}/{{sampleID}}/{{sampleID}}.metabat-bins")
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.metabat.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_metabat.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Metabat \n"

        sampleID=$(echo $(basename $(dirname {input.assembly})))
        echo -e "\nCreating temporary directory {root_path}/{config[path][scratch]}/metabat/${{sampleID}} ... "
        mkdir -p {root_path}/{config[path][scratch]}/metabat/${{sampleID}}

        echo -e "\nMoving into temporary directory {root_path}/{config[path][scratch]}/metabat/${{sampleID}} ... "
        cd {root_path}/{config[path][scratch]}/metabat/${{sampleID}}

        mkdir -p {output}

        cp {input.assembly} ./contigs.fasta.gz
        ls -l ./contigs.fasta.gz >> {log}

        cp {input.depth}/*all.depth ./
        ls -l ./*.all.depth >> {log}

        gunzip -f ./contigs.fasta.gz
        ls -l ./contigs.fasta >> {log}

        du -h ./contigs.fasta >> {log}
        du -h ./*.all.depth >> {log}

        /usr/bin/time -v metabat2 -i contigs.fasta \
            -a ${{sampleID}}.all.depth \
            -s {config[params][metabatMin]} \
            -v \
            --seed {config[params][seed]} \
            -t 0 \
            -m {config[params][minBin]} \
            -o ${{sampleID}}

        if compgen -G "*.fa" > /dev/null; then
            mv *.fa {output}
            echo "Moved .fa files to output directory" >> {log}
        else
            mkdir -p {output}
            echo "No .fa files found to move" >> {log}
        fi

        rm -rf {root_path}/{config[path][scratch]}/metabat/${{sampleID}}

        echo "Metabat done at $(date)."
        """
