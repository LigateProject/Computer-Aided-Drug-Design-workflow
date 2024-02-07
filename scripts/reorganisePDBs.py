from pathlib import Path as path
from argparse import ArgumentParser as cmd_line_parser
import sys

# read in file names from the command line
parser = cmd_line_parser(
    prog="reorganisePDBs.py",
    description="fuses modelled protein structure with crystal waters and ions and creates PDB files needed to calculate the RMSD between the experimental and the modelled structure",
)
parser.add_argument(
    "pdb1",
    type=path,
    help='Path to the original ".pdb" file from PDBbind 2020',
)
parser.add_argument(
    "pdb2",
    type=path,
    help='Path to the ".pdb" file containing the modelled structure',
)
parser.add_argument(
    "missingResidues",
    type=path,
    help='Path to the file "missingResidues.txt" containing a list of residues that had to be modelled',
)
parser.add_argument(
    "residuesNotRemodelled",
    type=path,
    help='Path to the file "residuesNotRemodelled.txt" containing a list of residues that promod3 refused to remodel',
)
args = parser.parse_args()

# read input files
try:
    f = open(args.pdb1, "r")
    lines1 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from PDB file "{args.pdb1}"')
    raise err

try:
    f = open(args.pdb2, "r")
    lines2 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from PDB file "{args.pdb2}"')
    raise err

try:
    f = open(args.missingResidues, "r")
    lines3 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from file "{args.missingResidues}"')
    raise err

try:
    f = open(args.residuesNotRemodelled, "r")
    lines4 = f.readlines()
    f.close()
except IOError as err:
    log_error(f'Unable to read from file "{args.residuesNotRemodelled}"')
    raise err


# determine position of missing residues
positionOfMissingResidues = dict()

for line in lines3:
    data = line.split()
    hit = -1
    while hit < len(data) - 7:
        hit = data.index("chain", hit+1)
        if data[hit + 1] not in positionOfMissingResidues:
            positionOfMissingResidues.update({data[hit + 1] : []})
        positionOfMissingResidues[data[hit + 1]].append(int(data[hit + 4].split(")")[0]))

# determine residues that were not remodelled
positionOfResiduesNotRemodelled = dict()

for line in lines4:
    data = line.split()
    hit = -1
    while hit < len(data) - 10:
        hit = data.index("chain", hit+1)
        if data[hit + 1] not in positionOfResiduesNotRemodelled:
            positionOfResiduesNotRemodelled.update({data[hit + 1] : []})
        positionOfResiduesNotRemodelled[data[hit + 1]].append([data[hit + 4] , data[hit + 7].split(")")[0]])


# read in file with modelled structure
modelledStructure = []

for line in lines2:
    if "END" not in line:
        modelledStructure.append(line)


# read in data from original PDBbind 2020 file
# accept ACE cap as residue
aminoAcids = ["ARG", "HIS", "LYS", "ASP", "GLU", "ASN", "GLN", "SER", "THR", "CYS", "GLY", "PRO", "ALA", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "VAL", "ACE"]

proteinAtomsNoHydrogens = []
waterAndIons = []

f = open("CONECT.txt", "w")
g = open("CONECT_old.txt", "w")

residueOffset = 0
blockedResidueNumber = "-1000000" # we cannot expect residue numbering to start from values > 0 => choose a negative number that is really unlikely to occur
ACECaps = dict()
for line in lines1:
    if "ATOM" in line or "HETATM" in line:
        if line[17:20] in aminoAcids:
            if line[17:20] == "ACE":
                if line[21] not in ACECaps:
                    ACECaps.update({line[21] : []})
                ACECaps[line[21]].append(line)
            if line[12] != "H" and line[13] != "H":
                # Exclude residues that have not been remodelled by promod3
                ## to keep the logic functional, we need an empty list as entry if there are no missing residues or residues that were not remodelled in the respective chain
                if line[21] not in positionOfMissingResidues:
                    positionOfMissingResidues.update({line[21] : []})
                if line[21] not in positionOfResiduesNotRemodelled:
                    positionOfResiduesNotRemodelled.update({line[21] : []})
                ## We are reading from the original PDB file here. A residue read in cannot belong to the missing residues.
                ## It may appear to belong there, though, if the residue numbering in the original PDB file is very unconventional.
                if [line[17:20], line[22:27].strip()] not in positionOfResiduesNotRemodelled[line[21]] and line[22:27] != blockedResidueNumber:
                    proteinAtomsNoHydrogens.append(line)
                    # ensure consecutive residue numbering starting from 1 (only start from 0 for ACE)
                    if (len(proteinAtomsNoHydrogens) == 1 or "TER" in proteinAtomsNoHydrogens[-1]) and (line[17:20] != "ACE"):
                        residueOffset += 1
                    elif proteinAtomsNoHydrogens[-1][22:27] != originalResidueNumber:
                        residueOffset += 1
                    originalResidueNumber = proteinAtomsNoHydrogens[-1][22:27]
                    if  "TER" not in proteinAtomsNoHydrogens[-1]:
                        if len(proteinAtomsNoHydrogens) > 1 and "TER" not in proteinAtomsNoHydrogens[-2]:
                            ## promod3 uses consecutive numbering to indicate gaps
                            if residueOffset in positionOfMissingResidues[proteinAtomsNoHydrogens[-1][21]]:
                                while residueOffset in positionOfMissingResidues[proteinAtomsNoHydrogens[-1][21]]:
                                    residueOffset += 1
                        proteinAtomsNoHydrogens[-1] = proteinAtomsNoHydrogens[-1][0:22] + (str(residueOffset)).rjust(4)  + " " + proteinAtomsNoHydrogens[-1][27:]
                ## however, promod3 uses the original numbering to indicate which residues it does not remodel
                ## TODO: fix for unconventional residue numbering including letters
                else:
                    if line[22:27] != blockedResidueNumber:
                        blockedResidueNumber = line[22:27]
                        residueOffset += 1
                        positionOfResiduesNotRemodelled[line[21]][positionOfResiduesNotRemodelled[line[21]].index([line[17:20], line[22:27].strip()])] = [line[17:20], (str(residueOffset)).rjust(4) + " "]
        else:
            waterAndIons.append(line)
    if "TER" in line:
        # there might be letters in the line after TER that mess with the logic of this script
        # add new line after "TER" to comply with PDB format
        proteinAtomsNoHydrogens.append(line[:3] + "\n")
        residueOffset = 0
        blockedResidueNumber = "-1000000"
    if "CONECT" in line:
        # the fixed format of PDB does not ensure the presence of spaces
        data = [line[0:6]]
        for i in range(6, len(line), 5):
            data.append(line[i:i+5])
        # update numbering of CONECT statement
        ## identify atoms belonging to that number (assuming that atom numbers are displayed completely)
        atomData = []
        atomDataOld = []
        for k in range(len(proteinAtomsNoHydrogens)):
            count = 0
            for j in range(len(data)):
                if (data[j] + " ") in proteinAtomsNoHydrogens[k][0:22]:
                    character = proteinAtomsNoHydrogens[k][proteinAtomsNoHydrogens[k][0:22].find(data[j] + " ") - 1]
                    if character == " " or character == "M":
                        count += 1
                        atomData.append([j])
                        atomDataOld.append([j])
            if count == 1:
                atomData[-1].append(proteinAtomsNoHydrogens[k][12:26])
                atomDataOld[-1].append(proteinAtomsNoHydrogens[k][12:26])
                atomDataOld[-1].append(proteinAtomsNoHydrogens[k][6:11])
            elif count > 1:
                atomData.pop()
                atomDataOld.pop()
                atomData.pop()
                atomDataOld.pop()
                if count > 2:
                    atomData.pop()
                    atomDataOld.pop()
        for k in range(len(waterAndIons)):
            count = 0
            for j in range(len(data)):
                if (data[j] + " ") in waterAndIons[k][0:22]:
                    character = waterAndIons[k][waterAndIons[k][0:22].find(data[j] + " ") - 1]
                    if character == " " or character == "M":
                        count += 1
                        atomData.append([j])
                        atomDataOld.append([j])
            if count == 1:
                atomData[-1].append(waterAndIons[k][12:26])
                atomDataOld[-1].append(waterAndIons[k][12:26])
                atomData[-1].append(waterAndIons[k][6:11])
                atomDataOld[-1].append(waterAndIons[k][6:11])
            elif count > 1:
                atomData.pop()
                atomDataOld.pop()
                atomData.pop()
                atomDataOld.pop()
                if count > 2:
                    atomData.pop()
                    atomDataOld.pop()
        ## in the PS dataset, input PDB files may contain CONECT statements involving H atoms
        ## they can be recongnised by the fact that the list of atoms found in proteinAtomsNoHydrogens
        ## is shorter than the original CONECT statement after subtracting the word CONECT and the new line
        if len(atomDataOld) < len(data) - 2:
            continue
        ## get new number
        for i in range(len(atomData)):
            if len(atomData[i]) == 2:
                for k in range(len(modelledStructure)):
                    if atomData[i][1] == modelledStructure[k][12:26]:
                        atomData[i].append(modelledStructure[k][6:11])
        ## manipulate string
        for i in range(len(atomData)):
            # Surprisingly many PDB files in the data set have sequence records that are not consistent with the actual coordinates
            # If we fail at finding atoms in the modelled structure, this mismatch is most likely the reason
            # Print warning here and examine those PDB files after running the script
            if len(atomData[i]) < 3:
                errorFile = open("atomsInConectStatementsNotFound", "w")
                errorFile.close()
                print("Can't find atoms listed in CONECT statements in the modelled structure. This is probably due to a script error.")
                sys.exit(0)
            data[atomData[i][0] + 1] = atomData[i][2]
        lineToBeAdded = "CONECT"
        for i in range(1, len(data)):
            lineToBeAdded += data[i].rjust(5)
        lineToBeAdded += "\n"
        waterAndIons.append(lineToBeAdded)
        ## keep track of CONECT statements for later checks
        ### disulfide bonds are encoded as CONECT statements with three atoms => remove C atom, only keep S atoms
        if len(atomData) == 3:
            indexToRemove = -1
            for k in range(len(atomData) - 1):
                for l in range(k + 1, len(atomData)):
                    difference = int(atomData[k][2]) - int(atomData[l][2])
                    if difference == -1:
                        indexToRemove = k
                    elif difference == 1:
                        indexToRemove = l
            if indexToRemove > -1:
                newAtomData = []
                for k in range(len(atomData)):
                    if k != indexToRemove:
                        newAtomData.append(atomData[k])
                atomData = newAtomData
        atomData.sort()
        for i in range(len(atomData)):
            f.write("%s %s" % (atomData[i][2], atomData[i][1]))
        f.write("\n")
        for i in range(1, len(atomData)):
            f.write("bond %s %s\n" % (atomData[0][2], atomData[i][2]))
        atomDataOld.sort()
        for i in range(len(atomDataOld)):
            g.write("%s %s" % (atomDataOld[i][2], atomDataOld[i][1]))
        g.write("\n")
        for i in range(1, len(atomDataOld)):
            g.write("bond %s %s\n" % (atomDataOld[0][2], atomDataOld[i][2]))

f.close()
g.close()

# Some PDB files contain OXT atoms before gaps. They break the logic of this script and have to be removed.
proteinAtomsNoHydrogensFiltered = []
for i in range(len(proteinAtomsNoHydrogens) - 1):
    if ("OXT" in proteinAtomsNoHydrogens[i]) and (not "TER" in proteinAtomsNoHydrogens[i + 1]) and (not "END" in proteinAtomsNoHydrogens[i + 1]):
        continue
    proteinAtomsNoHydrogensFiltered.append(proteinAtomsNoHydrogens[i])
proteinAtomsNoHydrogensFiltered.append(proteinAtomsNoHydrogens[-1])
proteinAtomsNoHydrogens = proteinAtomsNoHydrogensFiltered

# convert to integer for later use
for chain in positionOfResiduesNotRemodelled:
    for i in range(len(positionOfResiduesNotRemodelled[chain])):
        positionOfResiduesNotRemodelled[chain][i][1] = int(positionOfResiduesNotRemodelled[chain][i][1])

# write fused structure
f = open("fused.pdb", "w")

ACECapWritten = []

for i in range(len(modelledStructure)):
    # add back ACE caps with original coordinates as they cannot be modelled by promod3
    if modelledStructure[i][21] in ACECaps and int(modelledStructure[i][22:26]) == 1 and modelledStructure[i][21] not in ACECapWritten:
        for j in range(len(ACECaps[modelledStructure[i][21]])):
            f.write(ACECaps[modelledStructure[i][21]][j])
        ACECapWritten.append(modelledStructure[i][21])
    f.write(modelledStructure[i])
for i in range(len(waterAndIons)):
    f.write(waterAndIons[i])
f.write("END\n")

f.close()

# write PDB file of original PDB entry without hydrogens
f = open("original.pdb", "w")

for i in range(len(proteinAtomsNoHydrogens)):
    f.write(proteinAtomsNoHydrogens[i])
f.write("END\n")

f.close()

# write modelled structure without modelled parts
f = open("remodelled.pdb", "w")
g = open("missingAtoms.txt" , "w")

index = 0
increment = 0
previousChainOriginal = "1" # chains are labelled with letters; a number cannot be the name of the previous chain
previousChainModelled = "1"
residueWithNewAtomOrder1 = []
residueWithNewAtomOrder2 = []
residueWithNewAtomOrder3 = []
residueWithNewAtomOrder4 = []

for i in range(len(modelledStructure)):
    # we cannot assume that the order of chains is preserved
    ## we might need to start with the chain listed last in the modelled structure
    if (i + increment) >= len(modelledStructure):
        increment = 0
    if "ATOM" in modelledStructure[i + increment]:
        ## ACE caps are preserved in proteinAtomsNoHydrogens but not in modelledStructure
        while "ACE" in proteinAtomsNoHydrogens[index]:
            index += 1
        ## avoid out-of-range errors
        if "ATOM" in proteinAtomsNoHydrogens[index]:
            chainName = proteinAtomsNoHydrogens[index][21]
            ## only update the increment if we have a new chain
            ## (We only update index if the atom is found in the original structure
            ## such that chains are effectively equal in length even if residues were modelled)
            if chainName != previousChainOriginal:
                ## calculate increment
                if chainName != modelledStructure[i][21]:
                    for j in range(len(modelledStructure)):
                        if "ATOM" in modelledStructure[j]:
                            if chainName == modelledStructure[j][21]:
                                foundChainName = j
                                break
                    increment = foundChainName - i
                else:
                    if chainName == previousChainModelled:
                        for j in range(len(modelledStructure)):
                            if "ATOM" in modelledStructure[j]:
                                if chainName == modelledStructure[j][21]:
                                    foundChainName = j
                                    break
                        increment = foundChainName - i
                    else:
                        increment = 0
                ## update memory
                previousChainOriginal = chainName
                previousChainModelled = modelledStructure[i][21]
        # only print missing atoms if they do not belong to missing residues or residues that were not remodelled
        if int(modelledStructure[i + increment][22:26].strip()) not in positionOfMissingResidues[modelledStructure[i + increment][21]] and [modelledStructure[i + increment][17:20], int(modelledStructure[i + increment][22:26])] not in positionOfResiduesNotRemodelled[modelledStructure[i + increment][21]]:
            # we cannot assume that the ordering of atoms within a residue is preserved; however, atom names are unique within a residue
            if index < len(proteinAtomsNoHydrogens):
                ## we have a complete non-terminal residue
                if len(residueWithNewAtomOrder1) > 0:
                    if modelledStructure[i + increment][22:26] != residueWithNewAtomOrder4[-1]:
                        indicesUsed = []
                        for j in range(len(residueWithNewAtomOrder3)):
                            indicesUsed.append(residueWithNewAtomOrder2.index(residueWithNewAtomOrder3[j]))
                            f.write(residueWithNewAtomOrder1[indicesUsed[-1]])
                        for j in range(len(residueWithNewAtomOrder2)):
                            if j not in indicesUsed:
                                # I don't want to use a modelled structure just because the second terminal oxygen is missing
                                if residueWithNewAtomOrder1[j][13:16] != "OXT":
                                    g.write("%d " % int(residueWithNewAtomOrder1[j][6:11]))
                        residueWithNewAtomOrder1 = []
                        residueWithNewAtomOrder2 = []
                        residueWithNewAtomOrder3 = []
                        residueWithNewAtomOrder4 = []
                ## collect atoms belonging to the residue
                if (modelledStructure[i + increment][17:20] in proteinAtomsNoHydrogens[index] and modelledStructure[i + increment][21] in proteinAtomsNoHydrogens[index]) and modelledStructure[i + increment][22:26] in proteinAtomsNoHydrogens[index][22:26]:
                    residueWithNewAtomOrder1.append(modelledStructure[i + increment])
                    residueWithNewAtomOrder2.append(modelledStructure[i + increment][13:16].strip())
                    residueWithNewAtomOrder3.append(proteinAtomsNoHydrogens[index][13:16].strip())
                    residueWithNewAtomOrder4.append(modelledStructure[i + increment][22:26])
                    index += 1
                ### we may have more modelled atoms than atoms in the original crystal structure
                else:
                    if len(residueWithNewAtomOrder4) == 0:
                        errorFile = open("indexError", "w")
                        errorFile.close()
                        print("IndexError in reorganisePDBs.py. This is probably due to a subtle and hard-to-fix error during structure repair (gaps were most likely identified incorrectly).")
                        sys.exit(0)
                    if modelledStructure[i + increment][22:26] == residueWithNewAtomOrder4[-1]:
                        residueWithNewAtomOrder1.append(modelledStructure[i + increment])
                        residueWithNewAtomOrder2.append(modelledStructure[i + increment][13:16].strip())
                        residueWithNewAtomOrder4.append(modelledStructure[i + increment][22:26])
                ## we have a complete terminal residue
                if len(residueWithNewAtomOrder1) > 0:
                    # the modelled structure may have an additional OXT that is not listed in the last residue of the original PDB file
                    # it is safe to assume that this additional OXT is always listed at the end of the residue
                    if "TER" in modelledStructure[i + increment + 1] or ("OXT" in modelledStructure[i + increment + 1] and index == len(proteinAtomsNoHydrogens)):
                        indicesUsed = []
                        for j in range(len(residueWithNewAtomOrder3)):
                            indicesUsed.append(residueWithNewAtomOrder2.index(residueWithNewAtomOrder3[j]))
                            f.write(residueWithNewAtomOrder1[indicesUsed[-1]])
                        for j in range(len(residueWithNewAtomOrder2)):
                            if j not in indicesUsed:
                                # I don't want to use a modelled structure just because the second terminal oxygen is missing
                                if residueWithNewAtomOrder1[j][13:16] != "OXT":
                                    g.write("%d " % int(residueWithNewAtomOrder1[j][6:11]))
                        residueWithNewAtomOrder1 = []
                        residueWithNewAtomOrder2 = []
                        residueWithNewAtomOrder3 = []
                        residueWithNewAtomOrder4 = []
    elif ("TER" in modelledStructure[i + increment] and "TER" in proteinAtomsNoHydrogens[index]):
        f.write(modelledStructure[i + increment])
        index += 1
    # if the original and the modelled structure don't have the same order of atoms, we stick to the order of atoms in the original PDB file
    # to avoid IndexError, we leave the for loop when we have reached the end of the final residue of the original structure
    # and the modelled structure contains at most an OXT atom and a "TER" signal we have not yet looped over
    if i + increment + 2 < len(modelledStructure):
        if ("OXT" in modelledStructure[i + increment + 1] and "TER" in modelledStructure[i + increment + 2]) and index == len(proteinAtomsNoHydrogens):
            break
    elif i + increment + 1 < len(modelledStructure):
        if "TER" in modelledStructure[i + increment + 1] and index == len(proteinAtomsNoHydrogens):
            break

f.write("END\n")

f.close()
g.close()
