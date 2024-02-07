import argparse
import os

import ost
from ost import io
from ost import conop
from ost import seq
from ost.mol.alg import Molck, MolckSettings

from promod3 import modelling

import sys

def _parse_args():

    desc = ("Performs structure standardisation using the ProMod3 homology "
            "modelling toolbox") 

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument("structure_in", help = "Structure file - supported "
                        "formats: .pdb, .pdb.gz, .cif, .cif.gz")
    parser.add_argument("seqres", help = "SEQRES sequences in fasta format - "
                        "expect one sequence per structure chain that should "
                        "be normalized. Matching happens with sequence/chain "
                        "names. Structure chains without matching sequence "
                        "are dropped.")
    parser.add_argument("structure_out", help = "Cleaned structure output in "
                        "PDB format")
    return parser.parse_args()

def _load_structure(structure_path):
    """Read OST entity either from mmCIF or PDB."""

    if not os.path.exists(structure_path):
        raise Exception(f"file not found: {structure_path}")

    # Determine file format from suffix.
    ext = structure_path.split(".")
    if ext[-1] == "gz":
        ext = ext[:-1]
    if len(ext) <= 1:
        raise Exception(f"Could not determine format of file "
                        f"{structure_path}.")
    sformat = ext[-1].lower()

    # increase loglevel, as we would pollute the info log with weird stuff
    ost.PushVerbosityLevel(ost.LogLevel.Error)
    # Load the structure
    if sformat in ["mmcif", "cif"]:
        entity = io.LoadMMCIF(structure_path)
        if len(entity.residues) == 0:
            raise Exception(f"No residues found in file: {structure_path}")
    elif sformat in ["pdb"]:
        try:
            entity = io.LoadPDB(structure_path)
        except Exception as e:
            print("Could not normalise structure: an exception occurred while loading PDB file with Open structure.")
            h = open("normalisationErrorReadingPDB", "w")
            h.close()
            sys.exit(0)
        if len(entity.residues) == 0:
            raise Exception(f"No residues found in file: {structure_path}")
    else:
        raise Exception(f"Unknown/ unsupported file extension found for "
                        f"file {structure_path}.")
    # restore old loglevel and return
    ost.PopVerbosityLevel()
    return entity

def _load_seqres(seqres_path):
    """ Read sequence list in fasta format
    """
    if not os.path.exists(seqres_path):
        raise Exception(f"file not found: {seqres_path}")
    seqres = io.LoadSequenceList(seqres_path, format = 'fasta')
    return seqres

def _clean(ent):
    """ Structure cleanup

    Performs the following processing:

    * to be documented
    """
    ms = MolckSettings(rm_unk_atoms=True,
                       rm_non_std=True,
                       rm_hyd_atoms=True,
                       rm_oxt_atoms=True,
                       rm_zero_occ_atoms=False,
                       colored=False,
                       map_nonstd_res=True,
                       assign_elem=True)
    Molck(ent, conop.GetDefaultLib(), ms)
    
    # requirement of running the processor here is considered a bug and is not
    # necessary starting from OpenStructure 2.4.0
    processor = conop.RuleBasedProcessor(conop.GetDefaultLib())
    processor.Process(ent)
    return ent.Select("peptide=true")

def _normalize(ent, seqres):
    """ GOGOGO
    """

    # we really want to know that stuff...
    ost.PushVerbosityLevel(3)

    # make some noise on missing chains/sequences
    seqres_names = set([s.GetName() for s in seqres])
    chain_names = set([ch.GetName() for ch in ent.chains])
    for sn in seqres_names:
        if sn not in chain_names:
            ost.LogInfo(f"seqres with name {sn} without match in structure")
    for chn in chain_names:
        if chn not in seqres_names:
            ost.LogInfo(f"chain with name {chn} without match in seqres")

    aln_list = seq.AlignmentList()
    names = list()
    for s in seqres:
        name = s.GetName()
        if name in chain_names:
            chain_sel = ent.Select(f"cname={name}")
            try:
                aln = seq.alg.AlignToSEQRES(chain_sel, s, validate=False)
            except Exception as e:
                # enrich with some more info an re-raise
                print("Could not normalise structure due to alignment failure")
                h = open("normalisationError", "w")
                h.close()
                sys.exit(0)
            aln.AttachView(1, chain_sel)
            aln_list.append(aln)
            names.append(name)

    mhandle = modelling.BuildRawModel(aln_list, chain_names = names)
    model = modelling.BuildFromRawModel(mhandle)
    ost.PopVerbosityLevel()
    return model

def main():
    args = _parse_args()
    ent = _load_structure(args.structure_in)
    seqres = _load_seqres(args.seqres)
    ent = _clean(ent)
    ent = _normalize(ent, seqres)
    io.SavePDB(ent, args.structure_out)

if __name__ == '__main__':
    main()
