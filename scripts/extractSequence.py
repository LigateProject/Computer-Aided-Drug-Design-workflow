from pathlib import Path as path
from argparse import ArgumentParser as cmd_line_parser

from Bio.Align import PairwiseAligner
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Chain
pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)

import sys
import copy
import re
from typing import List, Dict

# read in file names from the command line
parser = cmd_line_parser(
    prog="extractSequence.py",
    description="extract and compare sequences from FASTA file",
)
parser.add_argument(
    "fasta1",
    type=path,
    help='Path to the self-made ".fasta" file',
)
parser.add_argument(
    "fasta2",
    type=path,
    help='Path to the official PDB ".fasta" file',
)
parser.add_argument(
    "pdbfile",
    type=path,
    help='Path to PDB file',
)
parser.add_argument(
    "pdbnohetatm",
    type=path,
    help='Path to PDB file without heteroatoms (needed to identify gaps in the structure)',
)
parser.add_argument(
    "gapsinfo",
    type=path,
    help='Path to file containing information about gaps in the structure',
)
args = parser.parse_args()

# read input files
try:
    f = open(args.fasta1, "r")
    lines1 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from self-made FASTA file "{args.fasta1}"')
    raise err

try:
    f = open(args.fasta2, "r")
    lines2 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from official PDB FASTA file "{args.fasta2}"')
    raise err

try:
    f = open(args.pdbfile, "r")
    lines3 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from PDB file "{args.pdbfile}"')
    raise err

try:
    f = open(args.pdbnohetatm, "r")
    f.close()
except IOError as err:
    log_error(f'Unable to read from PDB file with deleted hetatoms "{args.pdbnohetatm}"')
    raise err

try:
    f = open(args.gapsinfo, "r")
    f.close()
except IOError as err:
    log_error(f'Unable to retrieve information about gaps. No {args.gapsinfo} file')
    raise err

# Load gaps found by pdb-tools package and select only the ones based on structure. "Seq. Gap" identifies sequence gap and "Found" marks
# the line with final statistics.
struct_gaps = dict()
with open(args.gapsinfo, 'r') as gap_file:
    for line in gap_file:
        if "Seq. Gap" not in line and 'Found' not in line:
            struct_gap_line = line.strip().split()
            chain_gap = struct_gap_line[0][0]
            if chain_gap not in struct_gaps.keys():
                struct_gaps[chain_gap] = []
            start_gap = struct_gap_line[0][5:]
            end_gap = struct_gap_line[-1][5:]
            struct_gaps[chain_gap].append([int(start_gap), int(end_gap)])
            
# Divide sequence into fragments, which can be used to build a complete sequence with filled in gaps.
# "chain" variable should be a Bio.PDB.Chain.Chain() object. For more information please refer to https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
# "struct_gaps" variable is a dictionary containing gap coordinates (start and end residue numbers) for each chain.            
def get_seq_fragments_from_chain(chain: Chain.Chain, struct_gaps: Dict[str, List[int]]):
    aminoAcids = dict({"ARG" : "R", "HIS" : "H", "LYS" : "K", "ASP" : "D", "GLU" : "E", "ASN" : "N", "GLN" : "Q", "SER" : "S", "THR" : "T", "CYS" : "C", "GLY" : "G", "PRO" : "P", "ALA" : "A", "ILE" : "I", "LEU" : "L", "MET" : "M", "PHE" : "F", "TYR" : "Y", "TRP" : "W", "VAL" : "V"})
    res_list = [res.id[1] for res in chain]
    res_names = [aminoAcids[res.resname] for res in chain]

    seq_fragments = []
    # Adding fake gaps at the beginning and at the end to get sequence parts easily
    gap_list = [[None, res_list[0]]] + struct_gaps[chain.id] + [[res_list[-1], None]]
    check_segmentation = True
    for i in range(1, len(gap_list)):
        start_fragment = gap_list[i-1][1]
        end_fragment = gap_list[i][0]
        seq_fragment = ''
        for num, resid in enumerate(res_list):
            # The following if-clauses cover the cases when gaps appear in hypervariable regions with special numbering e.g. 100A, 100B, 100C.
            # Most of the tools parse such numbering simply as 100, 100, 100 etc.
            # We check if gap start and end residue numbers are the same. If true, we do not add residues with numbers equal to gap start/end.
            # This causes intended skipping of certain residues present in the structure, which should not happen with other types of sequences.
            # We set "check_segmentation" as False to prevent further error rising for this kind of sequences.
            # Case 1: both gaps have the same start and end coordinates, e.g. gap i-1 is [100, 100] and gap i is [102, 102]
            if gap_list[i-1][0] == gap_list[i-1][1] and  gap_list[i][0] == gap_list[i][1]:
                check_segmentation = False
                if resid > start_fragment and resid < end_fragment:
                    seq_fragment+=res_names[num]
            # Case 2: only gap i-1 has the same start and end coordinates, e.g. gap i-1 is [100, 100] and gap i is [102, 105]
            elif gap_list[i-1][0] == gap_list[i-1][1]:
                check_segmentation = False
                if resid > start_fragment and resid <= end_fragment:
                    seq_fragment+=res_names[num]
            # Case 3: only gap i has the same start and end coordinates, e.g. gap i-1 is [95, 100] and gap i is [102, 102]
            elif gap_list[i][0] == gap_list[i][1]:
                check_segmentation = False
                if resid >= start_fragment and resid < end_fragment:
                    seq_fragment+=res_names[num]
            # Case 4: both gaps have different start and end coordinates, e.g. gap i-1 is [95, 100] and gap i is [102, 105]
            else:
                if resid >= start_fragment and resid <= end_fragment:
                    seq_fragment+=res_names[num]
        seq_fragments.append(seq_fragment)
    return seq_fragments, check_segmentation


# read in self-made fasta file
chains1 = dict()
keys1 = []

for line in lines1:
    if ">" in line:
        keys1.append(line[1])
        chains1.update({keys1[-1] : ""})
    else:
        chains1[keys1[-1]] += line[0:-1]

# read in PDB fasta file (consider that several chains can be summarised in one line for multimers)
chains2 = dict()
keys2 = []

for line in lines2:
    if ">" in line:
        if "Chains" in line:
            keys2.append(line[line.find("Chains") + 7])
        elif "Chain" in line:
            keys2.append(line[line.find("Chain") + 6])
        chains2.update({keys2[-1] : ""})
    else:
        chains2[keys2[-1]] += line[0:-1]


# Naming conventions in FASTA files are very complicated
# => Identify matching chains on the basis of alignment scores
chains2New = dict()
for chain1 in chains1:
    highScore = -sys.maxsize
    bestMatch = "?"
    for chain2 in chains2:
        aligner = PairwiseAligner()
        aligner.gap_score = -10
        alignment = aligner.align(chains1[chain1], chains2[chain2])[0]
        score = alignment.score
        if score > highScore:
            highScore = score
            bestMatch = chain2
    chains2New.update({chain1 : chains2[bestMatch]})
chains2 = chains2New


# I do not know how to interpret "X" as sequence in an official FASTA file => error out gracefully
## Doing it here to only exit if the chains present in the crystal structure are affected
for chain in chains2:
    if "X" in chains2[chain]:
        print("Found X in FASTA sequence")
        h = open("xInFastaSequence", "w")
        h.close()
        sys.exit(0)

# Check if there are gaps in the protein structure. If there are no gaps, the chain is already complete.
if struct_gaps == {}:
    #check if the sequence is correct according to official fasta
    print('No structural gaps found in the structure')
    for  k in chains1.keys():
        check = (chains1[k] in chains2[k])
        if check == False:
                print("It is possible that there is a mismatch between the sequence of the PDB structure and the sequence in the official FASTA file, or chain names have not been assigned correctly")
                h = open("alignmentError", "w")
                h.close()
                sys.exit(0)
    completed_chains = chains1
else:
    # In case there are gaps we build the sequence inside them based on official fasta files
    print(f'Found gaps: {len(struct_gaps.keys())}')
    print(struct_gaps)
    model = pdb_parser.get_structure('model', args.pdbnohetatm)[0]
    completed_chains = {}
    for chain in model:
        #prepare a list of sequence fragments between gaps
        if chain.id in struct_gaps.keys():
            crystal, check_segmentation = get_seq_fragments_from_chain(chain, struct_gaps)
            #check if extraction is correct
            if check_segmentation and ''.join(crystal) != chains1[chain.id]:
                print('Error in PDB file processing. Sequence after segmentation is not equal to original sequence from pdb.')
                h = open("alignmentError", "w")
                h.close()
                sys.exit(0)
            fasta = chains2[chain.id]

            for j in range(len(crystal)-1):
                beforeGap = crystal[0]
                afterGap = crystal[j+1]

                i = 0

                while i <= len(fasta) - len(beforeGap) - len(afterGap):
                    searchString = beforeGap + i * "." + afterGap
                    temp = re.compile(searchString)
                    fullCrystal = temp.search(fasta)
                    i += 1
                    if fullCrystal is not None:
                        crystal[0] = fullCrystal.group(0)
                        break

            final_chain_fasta = crystal[0]
            
        else:
            final_chain_fasta = chains1[chain.id]

        if final_chain_fasta not in chains2[chain.id] or len(final_chain_fasta) < len(chains1[chain.id]):
            print('Could not fill gaps in the PDB structure with the correct sequence!')
            h = open("alignmentError", "w")
            h.close()
            sys.exit(0)
        completed_chains[chain.id] = final_chain_fasta

for k in chains1.keys():
    print(k)
    print('PDB seq     : ', chains1[k])
    print('Official seq: ', chains2[k])        
    print('Final seq   : ', completed_chains[k])

# write output file
with open("new.fasta", "w") as f_out:
    for i in completed_chains.keys():
        f_out.write(f'>{i}\n')
        sequence = completed_chains[i]
        [f_out.write("%s\n" % sequence[j:j+60]) for j in range(0, len(sequence), 60)]