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
        expand(f"{config['path']['root']}/{config['folder']['ecGEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX)


####---------------------####
#### 8. GEMs to ecGEMs
####---------------------####

rule ecGEM:
    input:
        f"{config['path']['root']}/{config['folder']['GEMs']}/{{sampleID}}"
    output:
        f"{config['path']['root']}/{config['folder']['ecGEMs']}/{{sampleID}}"
    shell:
        """
        # Create the scratch folder
        idvar=$(basename {output})
        echo $idvar
        mkdir -p {output}
        cd {output}

        # Extract taxonomic information
        cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][UniprotID]} ./
        cp {config[path][root]}/{config[folder][gtdbtk]}/$idvar/${{idvar}}_gtdbtk_summary.tsv ./
        cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][kegg]} ./
        
        #module load cray-R
        Rscript {config[path][root]}/{config[folder][scripts]}/{config[scripts][taxonomy]}
        
        # Only copy the files with existing taxonomic information
        awk 'NR>1 {{print $1}}' taxonomy.txt | while read bin; do
            mkdir -p $bin/
            
            for folder in 'code' 'data' 'models' 'output'; do 
                mkdir -p $bin/$folder;
            done

            cp $(dirname $(dirname {output}))/GEMs/$idvar/$bin.xml $bin/models/
            cp {config[path][root]}/{config[folder][proteinBins]}/$idvar/$bin.faa $bin/
            awk '/>/ {{name = $0; gsub(/>/,"",name); gsub(" #.*","",name); gsub (/[-.]/,"_",name) ; next}} {{print name"\t"$0}} ' $bin/$bin.faa > $bin/$bin.faa.txt;

            # copy the matlab files into each bin folder
            cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][makeECGEM]} ./
            cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][extractMet]} ./
            cp {config[path][root]}/{config[folder][scripts]}/{config[scripts][adaptertemplate]} $bin/
            

            #manage the adapter template file
            AdapterFileName=$(echo "$bin" | sed 's/\.//g')ecGECKOAdapter
            mv $bin/adapterTemplate.m $bin/$AdapterFileName.m

            protID=$(awk -F '\t' -v genome="$bin" '$0 ~ genome {{print $2}}' taxonomy.txt)
            echo "Uniprot ID:"
            echo $protID
            ncbiID=$(awk -F '\t' -v genome="$bin" '$0 ~ genome {{print $3}}' taxonomy.txt)
            echo "NCBI ID:"
            echo $ncbiID
            species=$(awk -F '\t' -v genome="$bin" '$0 ~ genome {{print $6}}' taxonomy.txt)
            echo "Species: "
            echo $species
            keggorg=$(awk -F '\t' -v genome="$bin" '$0 ~ genome {{print $4}}' taxonomy.txt)
            echo "KEGG org: "
            echo $keggorg

            currentpath=$(pwd)
            params_path=$(echo "fullfile('$currentpath', '$bin')")

            cd $bin
            pwd

            sed -i.bak "
                /classdef/s/classdef.*/classdef ${{AdapterFileName}} < ModelAdapter/;
                /function obj =/s/function obj.*/function obj = ${{AdapterFileName}}()/;
                /obj.params.path/s|obj.params.path.*|obj.params.path = fullfile('${{currentpath}}', '${{bin}}');|;
                /obj.params.uniprot.type/s/obj.params.uniprot.type.*/obj.params.uniprot.type = 'proteome';/;
                /obj.params.complex.taxonomicID/s/obj.params.complex.taxonomicID.*/obj.params.complex.taxonomicID = '${{ncbiID}}';/;
                /obj.params.uniprot.ID/s/obj.params.uniprot.ID.*/obj.params.uniprot.ID = '${{protID}}';/;
                /obj.params.convGEM/s|obj.params.convGEM.*|obj.params.convGEM = fullfile(obj.params.path, 'models', '${{bin}}.xml');|;
                /obj.params.org_name/s/obj.params.org_name.*/obj.params.org_name = '${{species}}';/;
                /obj.params.uniprot.reviewed/s/obj.params.uniprot.reviewed.*/obj.params.uniprot.reviewed = false;/;
                /obj.params.kegg.ID/s/obj.params.kegg.ID.*/obj.params.kegg.ID = '${{keggorg}}';/;
                /obj.params.enzyme_comp/s/obj.params.enzyme_comp.*/obj.params.enzyme_comp = 'cytosol';/
            " ${{AdapterFileName}}.m

            rm ${{AdapterFileName}}.m.bak

            cd ../

        done

        # Run make-gecko
        cd {output}
        ecGEMpath=$(echo $(dirname $(pwd))/$idvar)
        GEMpath=$(echo $(dirname $(dirname $ecGEMpath))/GEMs/$idvar)
        
        echo $GEMpath
        echo $ecGEMpath

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '${{GEMpath}}';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '${{ecGEMpath}}';|;
        " meta-gecko-extractMet.m

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '${{GEMpath}}';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '${{ecGEMpath}}';|;
        " meta-gecko.m

        rm meta-gecko.m.bak
        rm meta-gecko-extractMet.m.bak

        #/Applications/MATLAB_R2024b.app/toolbox/matlab/general/matlab.m -nodisplay -nosplash -nodesktop -logfile matlab.log/$idvar.txt < meta-gecko.m

        # # DLKcat needs to run locally

        # for bin in bin.1.p/; do
        #     dlkcat_result_nokcat = $(echo $bin/data/DLKcat.tsv)
        #     dlkcat_results_kcat = $(echo $bin/data/DLKcatOutput.tsv)
        #     dlkcat_py = {config[path][root]}/{config[folder][dlkcat]}/{config[scripts][dlkcat]}
        #     cd $(dirname $dlkcat_py)

        #     python3.10 $dlkcat_py $dlkcat_result_nokcat $dlkcat_results_kcat

        #     if [ -f $dlkcat_results_kcat ]; then
        #         rm $dlkcat_result_nokcat
        #         mv $dlkcat_results_kcat $dlkcat_result_nokcat
        #         echo 'DLKcat prediction completed'
        #     else   
        #         echo 'DLKcat encountered an error or it did not create any output file.'
        #     fi
        # done

        # /Applications/MATLAB_R2024b.app/toolbox/matlab/general/matlab.m -nodisplay -nosplash -nodesktop -logfile matlab_log2.txt < meta-gecko-part2.m
 

        echo -e "ecGECKO models have been generated for $idvar at $(date)\nSection ends\n"


        """
