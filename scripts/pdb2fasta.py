from pathlib import Path as path
from argparse import ArgumentParser as cmd_line_parser

import sys

# read in file names from the command line
parser = cmd_line_parser(
    prog="pdb2fasta.py",
    description="get FASTA sequence from PDB file",
)
parser.add_argument(
    "pdbfile",
    type=path,
    help='Path to the target ".pdb" file',
)
args = parser.parse_args()

try:
    f = open(args.pdbfile, "r")
    lines = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from PDB file "{args.pdbfile}"')
    raise err

# dictionary for one-letter code
aminoAcids = dict({"ARG" : "R", "HIS" : "H", "LYS" : "K", "ASP" : "D", "GLU" : "E", "ASN" : "N", "GLN" : "Q", "SER" : "S", "THR" : "T", "CYS" : "C", "GLY" : "G", "PRO" : "P", "ALA" : "A", "ILE" : "I", "LEU" : "L", "MET" : "M", "PHE" : "F", "TYR" : "Y", "TRP" : "W", "VAL" : "V", "HID" : "H", "HIE" : "H", "HIP" : "H", "LYN" : "K", "CYM" : "C", "CYX" : "C", "ASH" : "D","GLH" : "E"})

# convert PDB to FASTA files
sequence = dict()
previousResidue = "000"
countTER = 0
for line in lines:
    # no input PDB file contains the keyword "HETATM" in this data set
    # TODO: how to handle water molecules and ions in CADD workflow
    if "ATOM" in line:
        if line[21] not in sequence:
            sequence.update({line[21] : ""})
        if line[17:20] not in aminoAcids and line[17:20] != "ACE" and line[17:20] != "NME":
            h = open("nonStandardResidueFound", "w")
            h.close()
            print("Non-standard residue (%s) found in protein coordinates!" % (line[17:20]))
            sys.exit(0)
        if line[22:27] != previousResidue and line[17:20] != "ACE" and line[17:20] != "NME":
            sequence[line[21]] += aminoAcids[line[17:20]]
            previousResidue = line[22:27]
    if "TER" in line:
        countTER += 1

# write output file
#TODO: remove assumption about file naming
g = open(str(args.pdbfile).split("/")[-1].split(".")[0] + ".fasta", "w")
for chain in sequence:
    g.write("%s\n" % (">" + chain))
    [g.write("%s\n" % sequence[chain][j:j+60]) for j in range(0, len(sequence[chain]), 60)]
g.close()
