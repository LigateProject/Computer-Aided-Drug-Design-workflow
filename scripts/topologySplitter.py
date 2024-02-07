import sys

#print("Which file contains the topology to be split into two?")
filename = input()

inFile = open(filename, "r")
outFile1 = open("ffMOL.itp", "w")
outFile2 = open("ligandSingle.itp", "w")

lines = inFile.readlines()

check1 = False
check2 = False
check3 = False

for line in lines:
    if (len(line.split()) == 3):
        if (line.split()[1] == "atomtypes"):
            check1 = True
        elif (line.split()[1] == "moleculetype"):
            check1 = False
            check2 = True
    if check1:
        outFile1.write(line)
    elif check2:
        if check3:
            outFile2.write('MOL  3\n')
            check3 = False
        elif (line == lines[-2]):
            outFile2.write('#include "posre_ligand.itp"\n')
        else:
            outFile2.write(line)
        if (len(line.split()) == 3):
            if (line.split()[1] == "Name"):
                check3 = True

inFile.close()
outFile1.close()
outFile2.close()
