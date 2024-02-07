f = open("merged.itp", "r")
g = open("mergedConstantMass.itp", "w")

lines = f.readlines()

modify = False

for line in lines:
    data = line.split()
    if (modify and len(data) == 11):
        if (float(data[7]) < float(data[10])):
            mass = data[10]
        else:
            mass = data[7]
        g.write("%6s%12s%7s%7s%7s%7s%11s%11s%12s%11s%11s\n" % (data[0], data[1], data[2], data[3], data[4], data[5], data[6], mass, data[8], data[9], mass))
    else:
        g.write(line)
    if (len(data) > 1 and data[1] == "atoms"):
        modify = True
    elif (len(data) > 1 and data[1] == "bonds"):
        modify = False

f.close()
g.close()
