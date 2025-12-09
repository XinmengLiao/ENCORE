"""
ENCORE Downstream Analysis Pipeline
Integrated workflow: Binning -> Refinement -> Extraction -> Classification -> GEM Reconstruction
Sequence: CONCOCT -> MaxBin2 -> MetaBat2 -> Refinement -> Reassembly -> Extraction -> Taxonomy -> CarveMe -> GECKO3
"""

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

# ================================
# Include individual workflow modules in execution order
# ================================
# Binning modules
include: "concoct.smk"
include: "maxbin.smk"
include: "metabat.smk"
# Refinement and reassembly
include: "refinement.smk"
include: "reassembly.smk"
# Extraction and classification
include: "extraction.smk"
include: "taxonomy.smk"
# Metabolic modeling
include: "carveme.smk"
include: "GEMtoECGEM.smk"

# ================================
# Pipeline-level rule
# ================================
rule all:
    input:
        # Binning outputs
        expand(f"{root_path}/{config['folder']['concoct']}/{{sampleID}}/{{sampleID}}.concoct-bins", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['maxbin']}/{{sampleID}}/{{sampleID}}.maxbin-bins", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['metabat']}/{{sampleID}}/{{sampleID}}.metabat-bins", sampleID=SAMPLE_INDEX),
        # Refinement and reassembly
        expand(f"{root_path}/{config['folder']['refined']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['reassembled']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        # Extraction outputs
        expand(f"{root_path}/{config['folder']['dnaBins']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        expand(f"{root_path}/{config['folder']['proteinBins']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        # Taxonomy outputs
        expand(f"{root_path}/{config['folder']['classification']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        # GEM outputs
        expand(f"{root_path}/{config['folder']['GEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX),
        # ecGEM outputs
        expand(f"{root_path}/{config['folder']['ecGEMs']}/{{sampleID}}", sampleID=SAMPLE_INDEX)
    message:
        """
        This is the downstream analyses after 1. preprocessing.smk and 2. crossmapping.smk.
        Please ensure to complete the previous modules first, which could be run by:

        bash ENCORE.sh -o test/ -c test -d Toy_Dataset --preprocessing --crossmapping

        """