f = open("merged.itp", "r")
g = open("posre_ligand.itp", "w")

check = -2

lines = f.readlines()

g.write("[ position_restraints ]\n")
g.write("; atom  type      fx      fy      fz\n")

for line in lines:
    data = line.split()
    if check == 0:
        if len(data) > 1:
            if data[1] == "bonds":
                check = -2
            elif not 'h' in data[1] and not 'H' in data[1]:
                g.write("%6d%6d%6d%6d%6d\n" % (int(data[0]), 1, 1000, 1000, 1000))
    elif check == -1:
        check = 0
    elif len(data) > 1:
        if data[1] == "atoms":
            check = -1


f.close()
g.close()
