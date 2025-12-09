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
        expand(f"{root_path}/{config['folder']['GEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX)



####------------------------####
#### 6. GEMs Constructions
####------------------------####


#############
## Carveme ## (can run, but need to set other parameters first)
#############
rule carveme:
    input:
        bin = f"{root_path}/{config['folder']['proteinBins']}/{{sampleID}}"
    output:
        directory(f"{root_path}/{config['folder']['GEMs']}/{{sampleID}}")
    params:
        scratch_dir = f"{root_path}/{config['path']['scratch']}/{config['folder']['GEMs']}"
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.carveme.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_carveme.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Carveme \n"

        idvar=$(echo $(basename {output}))
        mkdir -p {output}
        mkdir -p {params.scratch_dir}/${{idvar}}
        cd {params.scratch_dir}/${{idvar}}

        cp {input.bin}/* ./
        medie_file=$(echo $(dirname {root_path})/{config[path][scripts]}/{config[DataforScripts][carveme]})
        cp $medie_file ./

        for name in *.faa; do
            binID=$(echo $name | sed 's/.faa//g')
            echo "Begin carving GEM for sample: ${{idvar}} protein bin: ${{binID}} ... "
            carve $name \
                -g {config[params][carveMedia]} \
                -v \
                --mediadb $medie_file \
                --fbc2 \
                -o ${{binID}}.xml 
        
            echo "${{idvar}} protein bin: ${{binID}} finished constructing GEMs at $(date)";
        done

        mkdir -p {output}
        mv *.xml {output}
        echo "Done carving GEM (Carveme) for all protein bins of Sample: ${{idvar}}. Now start to generate the statistic for GEMs."

        cd {output}
        while read model;do
            id=$(echo $(basename $model)|sed 's/.xml//g');
            mets=$(less $model| grep "species metaid="|cut -d ' ' -f 8|sort|uniq|wc -l);
            rxns=$(less $model|grep -c 'reaction metaid=');
            genes=$(less $model|grep -E 'fbc:geneProduct.*fbc:id='|wc -l);
            echo "Model: $id has $mets mets$, $rxns reactions, and $genes genes ... "
            echo "$id $mets $rxns $genes" >> GEMs.stats;
        done< <(find . -name "*.xml")

        mkdir -p {root_path}/{config[folder][stats]}/${{idvar}}
        mv GEMs.stats {root_path}/{config[folder][stats]}/${{idvar}}
        cd {root_path}/{config[folder][stats]}/${{idvar}}
        Rscript $(dirname {root_path})/{config[Rscripts][modelVis]}
        echo "GEMs statistic checking done at $(date)."
        """

