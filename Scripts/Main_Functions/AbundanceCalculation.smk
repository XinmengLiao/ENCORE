configfile: "ENCORE_config.yaml"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern(f"{config['path']['root']}/{config['folder']['data']}/*") 

rule all:
	input:
        	expand(f"{config['path']['root']}/{config['folder']['abundance']}/{{sampleID}}", sampleID=SAMPLE_INDEX)




rule abundance:
    input:
        bins = f"{config['path']['root']}/{config['folder']['reassembled']}/{{sampleID}}",
        R1 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{config['path']['root']}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
    output:
        directory(f"{config['path']['root']}/{config['folder']['abundance']}/{{sampleID}}")
    message: 
		"""
        Calculate bin abundance fraction using the following:

        binAbundanceFraction = ( X / Y / Z) * 1000000

        X = # of reads mapped to bin_i from sample_k
        Y = length of bin_i (bp)
        Z = # of reads mapped to all bins in sample_k

        Note: 1000000 scaling factor converts length in bp to Mbp
              Rule slightly modified for european datasets where input bins are in dna_bins_organized
              instead of metaWRAP reassembly output folder 

        
        Required tools: bwa, samtools
        """
	shell:
		"""
        mkdir -p {output}
        idvar=$(basename {input.bins})
        mkdir -p {config[path][scratch]}/${{idvar}}

        echo -e "\nCopying quality filtered paired end reads and generated MAGs to SCRATCHDIR ... "
        cd {config[path][scratch]}/${{idvar}}
        cp {input.R1} {input.R2} {input.bins}/${{idvar}}/reassembled_bins/* ./

        echo -e "\nConcatenating all bins into one FASTA file ... "
        cat *.fa > $(basename {input.bins}).fa

        echo -e "\nCreating bwa index for concatenated FASTA file ... "
        bwa index $(basename {input.bins}).fa

        echo -e "\nMapping quality filtered paired end reads to concatenated FASTA file with bwa mem ... "
        bwa mem -t {config[cores][abundance]} $(basename {input.bins}).fa \
            $(basename {input.R1}) $(basename {input.R2}) > $(basename {input.bins}).sam

        echo -e "\nConverting SAM to BAM with samtools view ... "
        samtools view -@ {config[cores][abundance]} -Sb $(basename {input.bins}).sam > $(basename {input.bins}).bam

        echo -e "\nSorting BAM file with samtools sort ... "
        samtools sort -@ {config[cores][abundance]} -o $(basename {input.bins}).sort.bam $(basename {input.bins}).bam

        echo -e "\nExtracting stats from sorted BAM file with samtools flagstat ... "
        samtools flagstat $(basename {input.bins}).sort.bam > map.stats

        echo -e "\nCopying sample_map.stats file to root/abundance/sample for bin concatenation and deleting temporary FASTA file ... "
        cp map.stats {output}/$(basename {input.bins})_map.stats
        rm $(basename {input.bins}).fa
        
        echo -e "\nRepeat procedure for each bin ... "
        for bin in *.fa;do

            echo -e "\nSetting up temporary sub-directory to map against bin $bin ... "
            mkdir -p $(echo "$bin"| sed "s/.fa//")
            mv $bin $(echo "$bin"| sed "s/.fa//")
            cd $(echo "$bin"| sed "s/.fa//")

            echo -e "\nCreating bwa index for bin $bin ... "
            bwa index $bin

            echo -e "\nMapping quality filtered paired end reads to bin $bin with bwa mem ... "
            bwa mem -t {config[cores][abundance]} $bin \
                ../$(basename {input.R1}) ../$(basename {input.R2}) > $(echo "$bin"|sed "s/.fa/.sam/")

            echo -e "\nConverting SAM to BAM with samtools view at $(date)... "
            samtools view -@ {config[cores][abundance]} -Sb $(echo "$bin"|sed "s/.fa/.sam/") > $(echo "$bin"|sed "s/.fa/.bam/")

            echo -e "\nSorting BAM file with samtools sort $(date)... "
            samtools sort -@ {config[cores][abundance]} -o $(echo "$bin"|sed "s/.fa/.sort.bam/") $(echo "$bin"|sed "s/.fa/.bam/")

            echo -e "\nExtracting stats from sorted BAM file with samtools flagstat at $(date) ... "
            samtools flagstat $(echo "$bin"|sed "s/.fa/.sort.bam/") > $(echo "$bin"|sed "s/.fa/.map/")

            echo -e "\nAppending bin length to bin.map stats file ... "
            echo -n "Bin Length = " >> $(echo "$bin"|sed "s/.fa/.map/")

            # Need to check if bins are original (megahit-assembled) or strict/permissive (metaspades-assembled)
            if [[ $bin == *.strict.fa ]] || [[ $bin == *.permissive.fa ]] || [[ $bin == *.s.fa ]] || [[ $bin == *.p.fa ]];then
                less $bin |grep ">"|cut -d '_' -f4|awk '{{sum+=$1}}END{{print sum}}' >> $(echo "$bin"|sed "s/.fa/.map/")
            else
                less $bin |grep ">"|cut -d '-' -f4|sed 's/len_//g'|awk '{{sum+=$1}}END{{print sum}}' >> $(echo "$bin"|sed "s/.fa/.map/")
            fi

            paste $(echo "$bin"|sed "s/.fa/.map/")

            echo -e "\nCalculating abundance for bin $bin ... "
            echo -n "$bin"|sed "s/.fa//" >> $(echo "$bin"|sed "s/.fa/.abund/")
            echo -n $'\t' >> $(echo "$bin"|sed "s/.fa/.abund/")

            X=$(less $(echo "$bin"|sed "s/.fa/.map/")|grep "mapped ("|awk -F' ' '{{print $1}}')
            Y=$(less $(echo "$bin"|sed "s/.fa/.map/")|tail -n 1|awk -F' ' '{{print $4}}')
            Z=$(less "../map.stats"|grep "mapped ("|awk -F' ' '{{print $1}}')
            awk -v x="$X" -v y="$Y" -v z="$Z" 'BEGIN{{print (x/y/z) * 1000000}}' >> $(echo "$bin"|sed "s/.fa/.abund/")
            
            paste $(echo "$bin"|sed "s/.fa/.abund/")
            
            echo -e "\nRemoving temporary files for bin $bin ... "
            rm $bin
            cp $(echo "$bin"|sed "s/.fa/.map/") {output}
            mv $(echo "$bin"|sed "s/.fa/.abund/") ../
            cd ..
            rm -r $(echo "$bin"| sed "s/.fa//")
        done

        echo -e "\nDone processing all bins, summarizing results into sample.abund file ... "
        cat *.abund > $(basename {input.bins}).abund

        echo -ne "\nSumming calculated abundances to obtain normalization value ... "
        norm=$(less $(basename {input.bins}).abund |awk '{{sum+=$2}}END{{print sum}}');
        echo $norm

        echo -e "\nGenerating column with abundances normalized between 0 and 1 ... "
        awk -v NORM="$norm" '{{printf $1"\t"$2"\t"$2/NORM"\\n"}}' $(basename {input.bins}).abund > abundance.txt

        rm $(basename {input.bins}).abund
        mv abundance.txt $(basename {input.bins}).abund

        mv $(basename {input.bins}).abund {output}
		"""
