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

# Paths of scripts and data in the config are relative to the config file folder (passed by ENCORE.sh)
CONFIG_DIR = config.get('CONFIG_DIR') or os.path.dirname(os.path.dirname(workflow.basedir))

def config_path(value):
    """Config value as a path: expand ~, resolve relative to the config folder ('' if empty)"""
    if not value:
        return ''
    value = os.path.expanduser(str(value).strip('"\''))
    return value if os.path.isabs(value) else os.path.normpath(os.path.join(CONFIG_DIR, value))

def matlab_str(value):
    return value.replace("'", "''")

# MATLAB commands run before meta-gecko.m: toolboxes from 'matlab:' in the config, then ENCORE's
# versions of loadConventionalGEM, makeEcModel, fuzzyKcatMatching, getECfromDatabase_* on top of
# the path (used instead of the GECKO originals), and the docker location for DLKcat
matlab_cfg = config.get('matlab') or {}
MATLAB_INIT = []
for key in ['raven', 'gecko']:
    if config_path(matlab_cfg.get(key)):
        MATLAB_INIT.append(f"addpath(genpath('{matlab_str(config_path(matlab_cfg.get(key)))}'));")
if config_path(matlab_cfg.get('gurobi')):
    MATLAB_INIT.append(f"addpath('{matlab_str(config_path(matlab_cfg.get('gurobi')))}');")
MATLAB_INIT.append(f"addpath('{matlab_str(os.path.join(os.path.dirname(workflow.basedir), 'InstallationScripts'))}');")
MATLAB_INIT.append("setenv('PATH', [getenv('PATH') ':/usr/local/bin:/opt/homebrew/bin:/Applications/Docker.app/Contents/Resources/bin']);")
MATLAB_INIT = "\n".join(MATLAB_INIT)

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
        uniprotid_file = config_path(config['DataforScripts']['UniprotID']),
        kegg_file = config_path(config['DataforScripts']['kegg']),
        taxonomy_script = config_path(config['Rscripts']['taxonomy']),
        proteinbins_path = f"{root_path}/{config['folder']['proteinBins']}",
        gems_path = f"{root_path}/{config['folder']['GEMs']}",
        classification_path = f"{root_path}/{config['folder']['classification']}",
        makeecgem_script = config_path(config['MatlabScripts']['makeECGEM']),
        extractmet_script = config_path(config['MatlabScripts']['extractMet']),
        adaptertemplate_script = config_path(config['MatlabScripts']['adaptertemplate']),
        matlab_init = MATLAB_INIT
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

        # Escape a value for the replacement part of a sed s||| command
        sed_escape() {{ printf '%s' "$1" | sed -e 's/[\\\\&|]/\\\\&/g'; }}

        # Extract taxonomic information
        cp "{params.uniprotid_file}" ./UniprotID.txt
        cp {params.classification_path}/${{idvar}}/${{idvar}}_gtdbtk_summary.tsv ./
        cp "{params.kegg_file}" ./kegglist.txt
        
        Rscript "{params.taxonomy_script}"

        # copy the matlab files into the sample folder
        cp "{params.makeecgem_script}" ./meta-gecko.m
        cp "{params.extractmet_script}" ./meta-gecko-extractMet.m
        
        # Only copy the files with existing taxonomic information
        awk 'NR>1 {{print $1}}' taxonomy.txt | while read bin; do
            if [[ ! -f {params.gems_path}/${{idvar}}/$bin.xml || ! -f {params.proteinbins_path}/${{idvar}}/$bin.faa ]]; then
                echo "$bin: GEM or protein file missing, skipped"
                continue
            fi

            for folder in 'code' 'data' 'models' 'output'; do 
                mkdir -p $bin/$folder;
            done

            cp {params.gems_path}/${{idvar}}/$bin.xml $bin/models/
            cp {params.proteinbins_path}/${{idvar}}/$bin.faa $bin/
            # gene<TAB>full sequence, one row per gene (the .faa is line-wrapped, so join all sequence lines of a protein)
            awk '/>/ {{if (name != "") print name"\\t"seq; name = $0; gsub(/>/,"",name); gsub(" #.*","",name); gsub (/[-.]/,"_",name); seq = ""; next}} {{gsub(/[ \\t\\r]/,"",$0); seq = seq $0}} END {{if (name != "") print name"\\t"seq}}' $bin/$bin.faa > $bin/$bin.faa.txt;

            #manage the adapter template file, class name = bin name without dots, e.g. bin.1.s -> bin1secGECKOAdapter
            AdapterFileName=$(echo "$bin" | sed 's/\\.//g')ecGECKOAdapter
            cp "{params.adaptertemplate_script}" $bin/$AdapterFileName.m

            # exact match on the bin column of taxonomy.txt
            protID=$(awk -F '\\t' -v genome="$bin" '$1 == genome {{print $2; exit}}' taxonomy.txt)
            echo "Uniprot ID: $protID"
            ncbiID=$(awk -F '\\t' -v genome="$bin" '$1 == genome {{print $3; exit}}' taxonomy.txt)
            echo "NCBI ID: $ncbiID"
            species=$(awk -F '\\t' -v genome="$bin" '$1 == genome {{print $6; exit}}' taxonomy.txt)
            echo "Species: $species"
            keggorg=$(awk -F '\\t' -v genome="$bin" '$1 == genome {{print $4; exit}}' taxonomy.txt)
            echo "KEGG org: $keggorg"

            currentpath=$(pwd)

            sed -i.bak "
                /classdef/s|classdef.*|classdef ${{AdapterFileName}} < ModelAdapter|;
                /function obj =/s|function obj.*|function obj = ${{AdapterFileName}}()|;
                /obj.params.path/s|obj.params.path.*|obj.params.path = fullfile('$(sed_escape "$currentpath")', '${{bin}}');|;
                /obj.params.uniprot.type/s|obj.params.uniprot.type.*|obj.params.uniprot.type = 'proteome';|;
                /obj.params.complex.taxonomicID/s|obj.params.complex.taxonomicID.*|obj.params.complex.taxonomicID = '$(sed_escape "$ncbiID")';|;
                /obj.params.uniprot.ID/s|obj.params.uniprot.ID.*|obj.params.uniprot.ID = '$(sed_escape "$protID")';|;
                /obj.params.convGEM/s|obj.params.convGEM.*|obj.params.convGEM = fullfile(obj.params.path, 'models', '${{bin}}.xml');|;
                /obj.params.org_name/s|obj.params.org_name.*|obj.params.org_name = '$(sed_escape "$species")';|;
                /obj.params.uniprot.reviewed/s|obj.params.uniprot.reviewed.*|obj.params.uniprot.reviewed = false;|;
                /obj.params.kegg.ID/s|obj.params.kegg.ID.*|obj.params.kegg.ID = '$(sed_escape "$keggorg")';|;
                /obj.params.enzyme_comp/s|obj.params.enzyme_comp.*|obj.params.enzyme_comp = 'cytosol';|
            " $bin/${{AdapterFileName}}.m

            rm $bin/${{AdapterFileName}}.m.bak

        done

        # Run make-gecko
        cd {output}
        ecGEMpath=$(pwd)
        GEMpath={params.gems_path}/$idvar
        sensitivityValue={config[params][FVAsensitivity]}
        
        echo $GEMpath
        echo $ecGEMpath

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '$(sed_escape "$GEMpath")';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '$(sed_escape "$ecGEMpath")';|;
        " meta-gecko-extractMet.m

        sed -i.bak "
            /GEMpath =/s|GEMpath =.*|GEMpath = '$(sed_escape "$GEMpath")';|;
            /ecGEMpath =/s|ecGEMpath =.*|ecGEMpath = '$(sed_escape "$ecGEMpath")';|;
            /sensitivity =/s|sensitivity =.*|sensitivity = ${{sensitivityValue}};|;
        " meta-gecko.m

        rm meta-gecko.m.bak
        rm meta-gecko-extractMet.m.bak

        # toolbox paths (config 'matlab:') and ENCORE's GECKO functions, then meta-gecko.m
        printf '%s\\n' "{params.matlab_init}" > encore_matlab_init.m
        cat encore_matlab_init.m meta-gecko.m | matlab -nodisplay -nosplash -nodesktop -logfile matlab.log

        echo -e "ecGECKO models have been generated for $idvar at $(date)\\nSection ends\\n"
        """
