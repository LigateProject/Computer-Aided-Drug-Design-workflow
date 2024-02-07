pathToForceField = 'amber99sb-ildn.ff'

f1 = open("topol_ligandInWater.top", "w")

f1.write('; Include forcefield parameters\n')
f1.write('#include "' + pathToForceField + '/forcefield.itp"\n')
f1.write('#include "ffMOL.itp"\n')
f1.write('\n')
f1.write('; Include ligand topology\n')
f1.write('#include "ligandSingle.itp"\n')
f1.write('\n')
f1.write('; Include water topology\n')
f1.write('#include "' + pathToForceField + '/tip3p.itp"\n')
f1.write('\n')
f1.write('; Include topology for ions\n')
f1.write('#include "' + pathToForceField + '/ions.itp"\n')
f1.write('\n')
f1.write('[ system ]\n')
f1.write('; Name\n')
f1.write('ligand in water\n')
f1.write('\n')
f1.write('[ molecules ]\n')
f1.write('; Compound        #mols\n')
f1.write('MOL 1\n')

f1.close()

f2 = open("../topol.top", "r")
f3 = open("topol_amber.top", "w")

lines = f2.readlines()

check = False
counter = 0
haveToPrint = True

for line in lines:
    if (line == "; Include forcefield parameters\n"):
        check = True
        counter = 0
    if check:
        if (line == "Protein\n"):
            f3.write("Protein-ligand complex in water\n")
            continue
        if (line == "; Include water topology\n"):
            f3.write('; Include ligand topology\n')
            f3.write('#include "ligandSingle.itp"\n')
            f3.write('\n')
        if "SOL" in line:
            f3.write('MOL 1\n')
            haveToPrint = False
            continue
        f3.write(line)
        counter += 1
        if (counter == 2):
            f3.write('#include "ffMOL.itp"\n')

if haveToPrint:
    f3.write('MOL 1\n')
else:
    for line in lines:
        if "SOL" in line:
            f3.write(line)

f2.close()
f3.close()
