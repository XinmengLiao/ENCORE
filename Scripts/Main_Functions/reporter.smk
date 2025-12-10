import os
import glob

configfile: os.path.join(workflow.basedir, "../../ENCORE_config.yaml")

# Build the root path - prefer path_root, then OUTPUT_DIR, then config value
root_path = config.get('path_root') or config.get('OUTPUT_DIR') or config['path']['root']

# Handle relative paths for root_path
if not os.path.isabs(root_path):
    root_path = os.path.join(os.path.dirname(workflow.basedir), root_path)


rule all:
    input:
        f"{root_path}/{config['folder']['reporter']}/reporter_metabolites.txt"


####--------------------------------------####
#### Reporter Metabolites Identification
####--------------------------------------####

rule ReporterMetabolite:
    input:
        network = f"{root_path}/{config['folder']['network']}/MMNetwork.txt",
        fva = f"{root_path}/{config['folder']['network']}/fva_all_nonzero.tsv",
        daa = config.get('DAA_FILE', '')
    output:
        f"{root_path}/{config['folder']['reporter']}/reporter_metabolites.txt"
    shell:
        """
        mkdir -p $(dirname {output})

        Rscript $(dirname {root_path})/{config[Rscripts][identifyRM]} {input.network} {input.daa} {input.fva} $(dirname {output})
        """