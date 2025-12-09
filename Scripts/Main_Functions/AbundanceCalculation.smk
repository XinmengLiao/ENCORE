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
        expand(f"{root_path}/{config['folder']['abundance']}/{{sampleID}}", sampleID=SAMPLE_INDEX)


rule abundance:
    input:
        bins = f"{root_path}/{config['folder']['reassembled']}/{{sampleID}}",
        R1 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R1_qfiltered.fastq.gz",
        R2 = f"{root_path}/{config['folder']['qfiltered']}/{{sampleID}}/{{sampleID}}_R2_qfiltered.fastq.gz"
    output:
        directory(f"{root_path}/{config['folder']['abundance']}/{{sampleID}}")
    params:
        scratch_dir = f"{root_path}/{config['path']['scratch']}"
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.abundance.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_abundance.log"
    shell:
        """
        mkdir -p {output}
        idvar=$(basename {input.bins})
        mkdir -p {params.scratch_dir}/abundance/${{idvar}}

        echo -e "\nCopying quality filtered paired end reads and generated MAGs to SCRATCHDIR ... "
        cd {params.scratch_dir}/abundance/${{idvar}}
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
            bin_name=$(echo "$bin"|sed "s/.fa/.bam/")
            samtools view -@ {config[cores][abundance]} -Sb $(echo "$bin"|sed "s/.fa/.sam/") > "$bin_name"

            echo -e "\nSorting BAM file with samtools sort $(date)... "
            sort_name=$(echo "$bin"|sed "s/.fa/.sort.bam/")
            samtools sort -@ {config[cores][abundance]} -o "$sort_name" "$bin_name"

            echo -e "\nExtracting stats from sorted BAM file with samtools flagstat at $(date) ... "
            map_file=$(echo "$bin"|sed "s/.fa/.map/")
            sort_bam=$(echo "$bin"|sed "s/.fa/.sort.bam/")
            samtools flagstat "$sort_bam" > "$map_file"

            echo -e "\nAppending bin length to bin.map stats file ... "
            echo -n "Bin Length = " >> "$map_file"

            if [[ $bin == *.strict.fa ]] || [[ $bin == *.permissive.fa ]] || [[ $bin == *.s.fa ]] || [[ $bin == *.p.fa ]];then
                less $bin |grep ">"| cut -d '_' -f4|awk '{{sum+=$1}}END{{print sum}}' >> "$map_file"
            else
                less $bin |grep ">"| cut -d '-' -f4|sed 's/len_//g'|awk '{{sum+=$1}}END{{print sum}}' >> "$map_file"
            fi

            paste "$map_file"

            echo -e "\nCalculating abundance for bin $bin ... "
            abund_file=$(echo "$bin"|sed "s/.fa/.abund/")
            echo -n "$bin"|sed "s/.fa//" >> "$abund_file"
            echo -n $'\t' >> "$abund_file"

            X=$(less "$map_file"|grep "mapped ("|awk -F' ' '{{print $1}}')
            Y=$(less "$map_file"|tail -n 1|awk -F' ' '{{print $4}}')
            Z=$(less "../map.stats"|grep "mapped ("|awk -F' ' '{{print $1}}')
            awk -v x="$X" -v y="$Y" -v z="$Z" 'BEGIN{{print (x/y/z) * 1000000}}' >> "$abund_file"

            paste "$abund_file"

            echo -e "\nRemoving temporary files for bin $bin ... "
            rm $bin
            cp "$map_file" {output}
            mv "$abund_file" ../
            cd ..
            bin_dir=$(echo "$bin"| sed "s/.fa//")
            rm -r "$bin_dir"
        done

        echo -e "\nDone processing all bins, summarizing results into sample.abund file ... "
        cat *.abund > $(basename {input.bins}).abund

        echo -ne "\nSumming calculated abundances to obtain normalization value ... "
        norm=$(awk '{{sum+=$2}} END {{print sum}}' $(basename {input.bins}).abund)
        echo $norm

        echo -e "\nGenerating column with abundances normalized between 0 and 1 ... "
        awk -v NORM="$norm" "{{printf '%s\t%s\t%f\n', $1, $2, $2/NORM}}" $(basename {input.bins}).abund > abundance.txt

        rm $(basename {input.bins}).abund
        mv abundance.txt $(basename {input.bins}).abund

        mv $(basename {input.bins}).abund {output}
        """