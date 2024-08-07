import sys

#print("Which file contains the Ligen poses?")
filename = input()
#print("Please specify the SMILES string of the ligand.")
smiles = input()
#print("Which file should contain the extracted pose?")
filename2 = input()
number = filename2.split("_")[1]
number = number.split("/")[0]

try:
    intNumber = int(number)
except ValueError:
    sys.exit("Pose number has to be an integer number.")

inFile = open(filename, "r")
outFile = open(filename2, "w")

# check for Ligen dummy particles first
forbiddenLines = []
forbiddenNumbers = []
forbiddenAtoms = 0
forbiddenBonds = 0
maxCount = 4
count = -1
subCount = 0
activation = 0
for line in inFile:
    # undo signalling via count if it isn't the correct ligand (SMILES string provided in line after <TRIPOS>MOLECULE)
    if len(line.split()) > 0:
        if line.split()[0] != smiles:
            if activation > 0:
                count = -1
        activation = 0
    # once we started reading a new pose, count the number of "TRIPOS" occurrences
    if subCount > 0:
        if line.find("<TRIPOS>") > -1:
            subCount += 1
            if line.find("COMMENT") > -1:
                maxCount += 1
    if (count == intNumber) and (subCount > 0) and (subCount < 4):
        lineToTest = line.split()
        # check atoms
        if "Du" in lineToTest:
            forbiddenNumbers.append(lineToTest[0])
            forbiddenLines.append(line)
            forbiddenAtoms += 1
        # check bonds
        if (len(lineToTest) > 2) and (lineToTest[1] in forbiddenNumbers or lineToTest[2] in forbiddenNumbers):
            forbiddenLines.append(line)
            forbiddenBonds += 1
    # a new pose is indicated by "<TRIPOS>MOLECULE"
    if line.find("<TRIPOS>MOLECULE") > -1:
        count += 1
        subCount = 1
        activation = 1
    # a pose ends with the fourth occurrence of "TRIPOS"
    if subCount == maxCount:
        subCount = 0

inFile.close()

inFile = open(filename, "r")

# let's hope babel doesn't require atom numbering in MOL2 files to be continuous
maxCount = 4
count = -1
subCount = 0
stringToRemember = ""
for line in inFile:
    # undo signalling via count if it isn't the correct ligand (SMILES string provided in line after <TRIPOS>MOLECULE)
    if len(line.split()) > 0:
        if line.split()[0] == smiles:
            if count == intNumber:
                outFile.write(stringToRemember)
        else:
            if activation > 0:
                count = -1
            stringToRemember = ""
        activation = 0
    # once we started reading a new pose, count the number of "TRIPOS" occurrences
    if subCount > 0:
        if line.find("<TRIPOS>") > -1:
            subCount += 1
            if line.find("COMMENT") > -1:
                maxCount += 1
    if (count == intNumber) and (subCount > 0) and (subCount < 4):
        if line not in forbiddenLines:
            # correct header (dummy particles are always listed at the end of the atom list)
            lineToTest = line.split()
            if lineToTest != []:
                if (lineToTest[0] in forbiddenNumbers) and (len(lineToTest) == 5):
                    outFile.write("%5d %5d %5d %5d %5d\n" % (int(lineToTest[0])-forbiddenAtoms, int(lineToTest[1])-forbiddenBonds, int(lineToTest[2]), int(lineToTest[3]), int(lineToTest[4])))
                    continue
            outFile.write(line)
    # a new pose is indicated by "<TRIPOS>MOLECULE"
    if line.find("<TRIPOS>MOLECULE") > -1:
        count += 1
        subCount = 1
        stringToRemember = line
        activation = 1
    # a pose ends with the fourth occurrence of "TRIPOS"
    if subCount == maxCount:
        subCount = 0

inFile.close()
outFile.close()
