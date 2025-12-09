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
        expand(f"{root_path}/{config['folder']['dnaBins']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['proteinBins']}/{{sampleID}}", sampleID=SAMPLE_INDEX)

################################
## protein and dna extraction ##
################################
rule extractBins:
    input:
        reassembledBins = directory(f"{root_path}/{config['folder']['reassembled']}/{{sampleID}}")
    output:
        dna_bin = directory(f"{root_path}/{config['folder']['dnaBins']}/{{sampleID}}"),
        protein_bin = directory(f"{root_path}/{config['folder']['proteinBins']}/{{sampleID}}")
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.extraction.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_extraction.log"
    shell:
        """
        echo "Now start to extract protein bins at $(date)."

        mkdir -p {output.protein_bin}
        idvar=$(echo $(basename {output.dna_bin}))
        cd {input.reassembledBins}/${{idvar}}
        rm -rf work_files

        for bin in reassembled_bins.checkm/bins/*;do
            var=$(echo $bin/genes.faa | sed 's|reassembled_bins/||g'|sed 's|reassembled_bins.checkm/bins/||'|sed 's|/genes||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
            cp $bin/*.faa {output.protein_bin}/$var;
        done

        echo -e "Protein Extraction done at $(date).\nNow will extract DNA bins."

        mkdir -p {output.dna_bin}

        echo -e "Begin copying and renaming dna fasta bins from reassembled_bins/ to dna_bins/ ...\n"
        cd {input.reassembledBins}/${{idvar}}

        echo "Copying bins from sample ${{idvar}} ... "
        for bin in reassembled_bins/*;do
            var=$(echo $bin| sed 's|reassembled_bins/||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
            cp $bin {output.dna_bin}/$var;
        done

        echo -e "DNA bins extraction finishes at $(date). Now will do gtdbtk taxonomic assignment.\n"
        """
