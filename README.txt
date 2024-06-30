master_final.sh can run all stages of the workflow in a fully automatic fashion. In case you would like to run workflow stages by hand, master_final.sh also documents the order in which scripts have to be run to satisfy all relevant dependencies.

PRACTICAL HINTS HOW TO RUN THE CADD WORKFLOW
- To run the CADD workflow, the following five environment variables have to be set: CADD_SOFTWARE_MODULES, CADD_SCRIPTS_DIR, CADD_INPUT_DIR, CADD_OUTPUT_DIR and CADD_LOCAL_DIR. For example:
  export CADD_SOFTWARE_MODULES=/home/software/modules # path to module files for all required software dependencies
  export CADD_SCRIPTS_DIR=/home/CADD/scripts # path to folder with the CADD workflow scripts
  export CADD_INPUT_DIR=/home/proteinTarget # path to folder with the input PDB file for the protein and the input SMILES strings or MOL2 files for the ligands
  export CADD_OUTPUT_DIR=/home/CADDV/run # path to folder to which the output of the CADD workflow should be written
  export CADD_LOCAL_DIR=/scratch/CADD # path to folder where the I/O intensive operations of the CADD workflow can be performed locally on the compute node, bypassing the file system

Note:
SLURM settings in this repository are specific to Leonardo and may need to be updated on other clusters.
