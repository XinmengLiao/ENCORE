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
        expand(f"{root_path}/{config['folder']['ecGEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX)


####---------------------####
#### 8. GEMs to ecGEMs
####---------------------####

rule ecGEM:
    input:
        f"{root_path}/{config['folder']['GEMs']}/{{sampleID}}"
    output:
        directory(f"{root_path}/{config['folder']['ecGEMs']}/{{sampleID}}")
    params:
        scratch_dir = f"{root_path}/{config['path']['scratch']}/{config['folder']['ecGEMs']}",
        uniprotid_file = os.path.join(os.path.dirname(workflow.basedir), config['DataforScripts'].get('UniprotID', 'UniprotID.txt')),
        kegg_file = os.path.join(os.path.dirname(workflow.basedir), config['DataforScripts'].get('kegg', 'kegglist.txt')),
        taxonomy_script = os.path.join(os.path.dirname(os.path.dirname(workflow.basedir)), config['Rscripts']['taxonomy']),
        proteinbins_path = f"{root_path}/{config['folder']['proteinBins']}",
        gems_path = f"{root_path}/{config['folder']['GEMs']}",
        classification_path = f"{root_path}/{config['folder']['classification']}",
        makeecgem_script = os.path.join(os.path.dirname(os.path.dirname(workflow.basedir)), config['MatlabScripts']['makeECGEM']),
        extractmet_script = os.path.join(os.path.dirname(os.path.dirname(workflow.basedir)), config['MatlabScripts']['extractMet']),
        adaptertemplate_script = os.path.join(os.path.dirname(os.path.dirname(workflow.basedir)), config['MatlabScripts']['adaptertemplate'])
    benchmark:
        f"{root_path}/{config['folder']['benchmarks']}/{{sampleID}}.ecgem.benchmark.txt"
    log:
        f"{root_path}/{config['folder']['logs']}/{{sampleID}}_ecgem.log"
    shell:
        """
        echo -e "$(date)\\nSection starts\\necGEM\\n"

        idvar=$(basename {output})
        mkdir -p {output}
        cd {output}

        # Extract taxonomic information
        cp {params.uniprotid_file} ./
        cp {params.classification_path}/${{idvar}}/${{idvar}}_gtdbtk_summary.tsv ./
        cp {params.kegg_file} ./
        
        Rscript {params.taxonomy_script}
        
        # Only copy the files with existing taxonomic information
        awk 'NR>1 {{print $1}}' taxonomy.txt | while read bin; do
            mkdir -p $bin/
            
            for folder in 'code' 'data' 'models' 'output'; do 
                mkdir -p $bin/$folder;
            done

            cp {params.gems_path}/${{idvar}}/$bin.xml $bin/models/
            cp {params.proteinbins_path}/${{idvar}}/$bin.faa $bin/
            awk '/>/ {{name = $0; gsub(/>/,"",name); gsub(" #.*","",name); gsub (/[-.]/,"_",name) ; next}} {{print name"\t"$0}} ' $bin/$bin.faa > $bin/$bin.faa.txt;

            # copy the matlab files into each bin folder
            cp {params.makeecgem_script} ./
            cp {params.extractmet_script} ./
            cp {params.adaptertemplate_script} $bin/

            #manage the adapter template file
            AdapterFileName=$(echo "$bin" | sed 's/\.//g')ecGECKOAdapter
            mv $bin/adapterTemplate.m $bin/$AdapterFileName.m

            protID=$(awk -F '\\t' -v genome="$bin" '$0 ~ genome {{print $2}}' taxonomy.txt)
            echo "Uniprot ID:"
            echo $protID
            ncbiID=$(awk -F '\\t' -v genome="$bin" '$0 ~ genome {{print $3}}' taxonomy.txt)
            echo "NCBI ID:"
            echo $ncbiID
            species=$(awk -F '\\t' -v genome="$bin" '$0 ~ genome {{print $6}}' taxonomy.txt)
            echo "Species: "
            echo $species
            keggorg=$(awk -F '\\t' -v genome="$bin" '$0 ~ genome {{print $4}}' taxonomy.txt)
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
        sensitivityValue={config[params][FVAsensitivity]}
        
        echo $GEMpath
        echo $ecGEMpath

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '${{GEMpath}}';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '${{ecGEMpath}}';|;
        " meta-gecko-extractMet.m

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '${{GEMpath}}';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '${{ecGEMpath}}';|;
            /sensitivity =/s|sensitivity =.*|sensitivity = ${{sensitivityValue}};|;
        " meta-gecko.m

        rm meta-gecko.m.bak
        rm meta-gecko-extractMet.m.bak

        matlab -nodisplay -nosplash -nodesktop -logfile matlab.log/$idvar.txt < meta-gecko.m

        echo -e "ecGECKO models have been generated for $idvar at $(date)\\nSection ends\\n"
        """
