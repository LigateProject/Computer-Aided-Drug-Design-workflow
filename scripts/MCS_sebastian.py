#!/bin/env python3
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS, AllChem
import sys
from csv import DictWriter


def find_similaries_from_smiles(smiles: list):
    return find_similarities_from_mol(
        [Chem.MolFromSmiles(mol_smi) for mol_smi in smiles]
    )


def find_similarities_from_mol(molecule_list: list):
    # add the index of the molecule as input property
    for index, molecule in enumerate(molecule_list):
        molecule.SetProp("SimilarityIndex", str(index))

    # compute the similarity matrix (is actually a triangular matrix)
    MCS_pair_list = []
    for index, mol1 in enumerate(molecule_list):
        for mol2 in molecule_list[index + 1 :]:
            # compute the similarity scores
            similarity_RDK = DataStructs.TanimotoSimilarity(
                Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2)
            )
            similarity_RDK = round(similarity_RDK, 3)
            similarity_MORGAN = DataStructs.TanimotoSimilarity(
                AllChem.GetMorganFingerprint(mol1,2), AllChem.GetMorganFingerprint(mol2,2))
            similarity_MORGAN=round(similarity_MORGAN, 3)
            
            # compute max common subgraph on the molecules
            fmcs_res = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=True)
            mol_fmcs_res = Chem.MolFromSmarts(fmcs_res.smartsString)

            # append this information to the list
            MCS_pair_list.append(
                {
                    "MolIndex1": int(mol1.GetProp("SimilarityIndex")),
                    "MolIndex2": int(mol2.GetProp("SimilarityIndex")),
                    "TanimotoSimilarityRdk": similarity_RDK,
                    "TanimotoSimilarityMorgan": similarity_MORGAN,
                    "NumAtomsMol1": mol1.GetNumAtoms(),
                    "NumAtomsMol2": mol2.GetNumAtoms(),
                    "NumCommonAtoms": fmcs_res.numAtoms,
                    "NumCommonBonds": fmcs_res.numBonds,
                    "CommonAtomsMol1": mol1.GetSubstructMatch(mol_fmcs_res),
                    "CommonAtomsMol2": mol2.GetSubstructMatch(mol_fmcs_res),
                }
            )
    return MCS_pair_list


if __name__ == "__main__":
    # read the SMILES from the standard input
    molecule_list = [smiles.rstrip() for smiles in sys.stdin if smiles]

    # perform the elaboration
    similarity_matrix = find_similaries_from_smiles(molecule_list)

    # perform a small post-process to pretty-print the result
    for entry in similarity_matrix:
        entry["CommonAtomsMol1_str"] = "|".join(
            [str(x) for x in entry["CommonAtomsMol1"]]
        )
        entry["CommonAtomsMol2_str"] = "|".join(
            [str(x) for x in entry["CommonAtomsMol2"]]
        )
        entry["MolSMILES1"] = molecule_list[entry["MolIndex1"]]
        entry["MolSMILES2"] = molecule_list[entry["MolIndex2"]]
        del entry["MolIndex1"]
        del entry["MolIndex2"]
        del entry["CommonAtomsMol1"]
        del entry["CommonAtomsMol2"]

    # write the output data in CSV
    fields = [
        "MolSMILES1",
        "MolSMILES2",
        "TanimotoSimilarityRdk",
        "TanimotoSimilarityMorgan",
        "NumAtomsMol1",
        "NumAtomsMol2",
        "NumCommonAtoms",
        "NumCommonBonds",
        "CommonAtomsMol1_str",
        "CommonAtomsMol2_str",
    ]
    writer = DictWriter(sys.stdout, fieldnames=fields, delimiter=" ")
    writer.writeheader()
    for row in similarity_matrix:
        writer.writerow(row)
