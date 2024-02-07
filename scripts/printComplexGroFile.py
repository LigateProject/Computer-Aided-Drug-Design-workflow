class ComplexGroFilePrinter:

    def __init__(self, filenameA, filenameB):
        f1 = open(filenameA, "r")
        f2 = open(filenameB, "r")

        # read in .gro files
        self.proteinGro = []
        self.mergedGro = []
        groFiles = [f1, f2]
        groFilesInMemory = [self.proteinGro, self.mergedGro]
        index = 0
        index2 = 3
        indices = [[2, 3, 6], [3, 4, 7], [2, 2, 5]]
        for f in groFiles:
            groFilesInMemory[index].append("Protein in complex with merged ligand\n")
            lines = f.readlines()
            counter = 0
            for line in lines:
                if (counter > 0):
                    groFilesInMemory[index].append(line.split())
                    if (len(groFilesInMemory[index][-1]) == 1):
                        groFilesInMemory[index][-1] = int(groFilesInMemory[index][-1][0])
                    elif (line == lines[-1]):
                        groFilesInMemory[index][-1] = [float(groFilesInMemory[index][-1][i]) for i in range(len(groFilesInMemory[index][-1]))]
                    else:
                        if (len(groFilesInMemory[index][-1]) == 6):
                            index2 = 0
                        elif (len(groFilesInMemory[index][-1]) == 7):
                            index2 = 1
                        elif (len(groFilesInMemory[index][-1]) == 5):
                            index2 = 2
                        listToAppend = []
                        for i in range(indices[index2][0]):
                            listToAppend.append(groFilesInMemory[index][-1][i])
                        if (index2 != 2):
                            listToAppend.append(int(groFilesInMemory[index][-1][indices[index2][0]]))
                        for i in range(indices[index2][1], indices[index2][2]):
                            listToAppend.append(float(groFilesInMemory[index][-1][i]))
                        groFilesInMemory[index][-1] = listToAppend
                counter += 1
            index += 1

        f1.close()
        f2.close()

    def listToStringConverter(self, inputList):
        string = ""
        if (len(inputList) == 6):
            indices = [2, 3, 6]
            string += inputList[0].rjust(8)
            string += inputList[1].rjust(7)
            string += "{:5d}".format(inputList[indices[0]])
        elif (len(inputList) == 7):
            indices = [3, 4, 7]
            for i in range(3):
                string += inputList[i].rjust(5)
            string += "{:5d}".format(inputList[indices[0]])
        elif (len(inputList) == 5):
            indices = [2, 2, 5]
            string += inputList[0].rjust(8)
            string += inputList[1].rjust(12)
        for i in range(indices[1], indices[2]):
            string += "{:8.3f}".format(inputList[i])
        string += "\n"
        return string

    def printComplexGroFile(self, filename):
        f1 = open(filename, "w")

        f1.write(self.proteinGro[0])
        f1.write("{:5d}\n".format(self.proteinGro[1] + self.mergedGro[1]))
        index = 0
        for i in range(2, len(self.proteinGro) - 1):
            # stop before printing solvent
            if ("SOL" in self.proteinGro[i][0]) or ("HOH" in self.proteinGro[i][0]):
                index = i
                break
            f1.write(self.listToStringConverter(self.proteinGro[i]))
        for j in range(2, len(self.mergedGro) - 1):
            if (len(self.mergedGro[j]) == 6):
                self.mergedGro[j][2] += self.proteinGro[1]
            elif (len(self.mergedGro[j]) == 7):
                self.mergedGro[j][3] += self.proteinGro[1]
            f1.write(self.listToStringConverter(self.mergedGro[j]))
        # print solvent from protein PDB:
        if index != 0:
            for i in range(index, len(self.proteinGro) - 1):
                f1.write(self.listToStringConverter(self.proteinGro[i]))
        string = ""
        for i in range(len(self.proteinGro[-1])):
            string += "{:11.5f}".format(self.proteinGro[-1][i])
            if i < len(self.proteinGro[-1]) - 1:
                string += " "
        string += "\n"
        f1.write(string)

        f1.close()

#print("Which gro file contains the structure of the protein?")
filenameP = input()
#print("Which gro file contains the structure of the hybrid ligand?")
filenameHL = input()
#print("Which file should contain the structure of the complex?")
filenameM = input()

cgfp = ComplexGroFilePrinter(filenameP, filenameHL)
cgfp.printComplexGroFile(filenameM)
