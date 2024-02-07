import fileinput
import sys

#print("Which .gro file should be modified?")
filename = input()
#print("Desired change in box size?")
shift = float(input())

groFile = open(filename, "r")

lines = groFile.readlines()
stringToBeReplaced = lines[-1]
stringToBeReplacedList = stringToBeReplaced.split()
stringToReplace = ""
for i in range(len(stringToBeReplacedList)):
    stringToReplace += "{:11.5f}".format(float(stringToBeReplacedList[i]) + shift)
    if (i < len(stringToBeReplacedList) - 1):
        stringToReplace += " "
stringToReplace += "\n"
groFile.close()

for line in fileinput.input(filename, inplace = 1):
    if stringToBeReplaced in line:
        line = line.replace(stringToBeReplaced, stringToReplace)
    sys.stdout.write(line)
