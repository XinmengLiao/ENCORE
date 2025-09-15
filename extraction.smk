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

SAMPLE_INDEX = get_sampleids_from_path_pattern('Study_dataset/*')


rule all:
	input:
        	expand("{home}/dna_bins/{sampleID}",sampleID=SAMPLE_INDEX, home=HOME), ## Bin Refine
        	expand("{home}/protein_bins/{sampleID}",sampleID=SAMPLE_INDEX, home=HOME) ## Bin Reassemble

################################
## protein and dna extraction ##
################################
rule extractBins:
    input:reassembledBins = directory('{home}/reassembled_bins/{sampleID}')
    output:
        dna_bin = directory('{home}/dna_bins/{sampleID}'),
        protein_bin = directory('{home}/protein_bins/{sampleID}')
    envmodules:
        "bioinfo-tools",
        "conda"
    benchmark:
        '{home}/benchmark/{sampleID}.binReassemble.benchmark.txt'
    log:
        "{home}/logs/{sampleID}/{sampleID}_binreassemble.log"
    shell:
        """      
        echo "Now start to extrat protein bins at $(date). "

        mkdir -p {output.protein_bin}
	    idvar=$(echo $(basename {output.dna_bin}))
        cd {input.reassembledBins}/${{idvar}}
        rm -rf work_files

        for bin in reassembled_bins.checkm/bins/*;do
           var=$(echo $bin/genes.faa | sed 's|reassembled_bins/||g'|sed 's|reassembled_bins.checkm/bins/||'|sed 's|/genes||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
           cp $bin/*.faa {output.protein_bin}/$var;
        done

        echo -e "Protein Extraction done at $(date).\nNow will extract DNA bins. "

        # Make sure dnaBins folder exists
        mkdir -p {output.dna_bin}

        # Copy files
        echo -e "Begin copying and renaming dna fasta bins from reassembled_bins/ to dna_bins/ ... \n"
        cd {input.reassembledBins}/${{idvar}}

        # extract DNA bins
        echo "Copying bins from sample ${{idvar}} ... "
        for bin in reassembled_bins/*;do
            # Loop through each bin
            var=$(echo $bin| sed 's|reassembled_bins/||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
            cp $bin {output.dna_bin}/$var;
        done

        echo -e "DNA bins extraction finishes a $(date). Now will do gtdbtk taxonomic assignment. \n" 

        """
