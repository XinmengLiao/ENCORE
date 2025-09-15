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

SAMPLE_INDEX = get_sampleids_from_path_pattern('dataset2/*')


rule all:
	input:
        	expand("{home}/concoct/{sampleID}/cov", home = HOME, sampleID = SAMPLE_INDEX), ## CrossMap
        	expand("{home}/metabat/{sampleID}/cov",home = HOME, sampleID = SAMPLE_INDEX),
        	expand("{home}/maxbin/{sampleID}/cov", home = HOME, sampleID = SAMPLE_INDEX)


####--------------------------------------####
#### 3. CrossMapSeries (no problem, run smoothly)
####--------------------------------------####
# The final outputs are coverage tables 
# In this step, the bwa is processed 


rule crossMapSeries:
    input:
        contigs = "{home}/megahit_assembled/{sampleID}/contigs.fasta.gz"
    output:
        concoct = directory("{home}/concoct/{sampleID}/cov"),
        metabat = directory("{home}/metabat/{sampleID}/cov"),
        maxbin = directory("{home}/maxbin/{sampleID}/cov")
    benchmark:
        "{home}/benchmark/{sampleID}.crossMapSeries.benchmark.txt"
    threads:128
    envmodules:
        "bioinfo-tools"
    shell:
        """
        echo -e "$(date)\nSection starts\n CrossMap \n"

        # Actiavte the conda environment 
        #set +e
        #set +u 
        #source activate metagem 
        #set -u 
        #set -e

        # Create output folders
        mkdir -p {output.concoct}
        mkdir -p {output.metabat}
        mkdir -p {output.maxbin}

	# Make job specific scratch dir
        idvar=$(echo $(basename $(dirname {output.concoct})))
        echo -e "\nCreating temporary directory {config[path][root]}/{config[folder][crossMap]}/${{idvar}} ... "
        #mkdir -p /cfs/klemming/scratch/x/xinmengl/crossMap/${{idvar}}
	mkdir -p {config[path][root]}/{config[folder][crossMap]}/${{idvar}}

        # Move into scratch dir
        cd {config[path][root]}/{config[folder][crossMap]}/${{idvar}}
	#cd /cfs/klemming/scratch/x/xinmengl/crossMap/${{idvar}}

        # Copy files
        cp {input.contigs} ./

        # Define the focal sample ID, fsample: 
        # The one sample's assembly that all other samples' read will be mapped against in a for loop
        fsampleID=$(echo $(basename $(dirname {input.contigs})))
        echo -e "\nFocal sample: ${{fsampleID}} ... "

        echo "Renaming and unzipping assembly ... "
        mv $(basename {input.contigs}) $(echo ${{fsampleID}}|sed 's/$/.fa.gz/g')
        gunzip -f $(echo ${{fsampleID}}|sed 's/$/.fa.gz/g')

        echo -e "\nIndexing assembly ... "
        bwa index ${{fsampleID}}.fa

	module load samtools
        
        for folder in {config[path][root]}/{config[folder][qfiltered]}/*;do 

                id=$(basename $folder)

                echo -e "\nCopying sample ${{id}} to be mapped against the focal sample ${{fsampleID}} ..."
                mkdir -p {config[path][root]}/{config[folder][crossMap]}/${{idvar}}/
		cp $folder/*.gz {config[path][root]}/{config[folder][crossMap]}/${{idvar}}/

                echo -e "\nMapping sample to assembly ... "
                rm -rf ${{id}}.sam  #make sure the output file not exist
                bwa mem -t 128 ${{fsampleID}}.fa *.fastq.gz > ${{id}}.sam
                
                echo -e "\nConverting SAM to BAM with samtools view ... " 
                rm -rf ${{id}}.bam
                samtools view -@ 128 -Sb ${{id}}.sam > ${{id}}.bam

                echo -e "\nSorting BAM file with samtools sort ... " 
                
                samtools sort -@ 128 -o ${{id}}.sort ${{id}}.bam

                echo -e "\nRunning jgi_summarize_bam_contig_depths script to generate contig abundance/depth file for maxbin2 input ... "
                jgi_summarize_bam_contig_depths --outputDepth ${{id}}.depth ${{id}}.sort

                echo -e "\nMoving depth file to sample ${{fsampleID}} maxbin2 folder ... "
                mv ${{id}}.depth {output.maxbin}

                echo -e "\nIndexing sorted BAM file with samtools index for CONCOCT input table generation ... " 
                samtools index ${{id}}.sort
		
		rm *.fastq.gz
		rm ${{id}}.sam ${{id}}.bam
        done
        
        nSamples=$(ls {config[path][root]}/{config[folder][qfiltered]}|wc -l)
        echo -e "\nDone mapping focal sample ${{fsampleID}} agains ${{nSamples}} samples in dataset folder."

        echo -e "\nRunning jgi_summarize_bam_contig_depths for all sorted bam files to generate metabat2 input ... "
        jgi_summarize_bam_contig_depths --outputDepth ${{fsampleID}}.all.depth *.sort

        echo -e "\nMoving input file ${{fsampleID}}.all.depth to ${{fsampleID}} metabat2 folder... "
        mv ${{fsampleID}}.all.depth {output.metabat}

        echo -e "Done. \nCutting up contigs to 10kbp chunks, not to be used for mapping!"
        cut_up_fasta.py -c 10000 -o 0 -m ${{fsampleID}}.fa -b assembly_c10k.bed > assembly_c10k.fa

        echo -e "\nSummarizing sorted and indexed BAM files with concoct_coverage_table.py to generate CONCOCT input table ... " 
        concoct_coverage_table.py assembly_c10k.bed *.sort > coverage_table.tsv

        echo -e "\nMoving CONCOCT input table to ${{fsampleID}} concoct folder"
        mv coverage_table.tsv {output.concoct}

        echo -e "\nRemoving intermediate sorted bam files ... "
        rm *.sort*

        echo "CrossMap done at $(date). "
        """


