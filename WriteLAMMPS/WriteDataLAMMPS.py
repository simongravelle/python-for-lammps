def write(data):   
    f = open("data.lammps", "w")
    f.write('# LAMMPS data file \n\n')
    f.write(str(data.Natoms)+' atoms\n')
    f.write(str(data.Nbonds)+' bonds\n')
    f.write(str(data.Nangles)+' angles\n')
    f.write(str(data.Ndihedrals)+' dihedrals\n')
    f.write(str(data.Nimpropers)+' impropers\n')
    f.write('\n')
    f.write(str(data.Tatoms)+' atom types\n')
    f.write(str(data.Tbonds)+' bond types\n')
    f.write(str(data.Tangles)+' angle types\n')
    f.write(str(data.Tdihedrals)+' dihedral types\n')
    f.write(str(data.Timpropers)+' improper types\n')
    f.write('\n')
    f.write('-25 25 xlo xhi\n')
    f.write('-25 25 ylo yhi\n')
    f.write('-25 25 zlo zhi\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for nlin in range(len(data.atoms.T)):
        newline = data.atoms.T[nlin]
        for col in range(len(newline)):
            if col < 3:
                f.write(str(int(newline[col]))+' ')
            else :
                f.write(str(newline[col])+' ')
        f.write('\n')
    f.write('\n')
    f.write('Bonds\n')
    f.write('\n')
    for nlin in range(len(data.bonds.T)):
        newline = data.bonds.T[nlin]
        for col in range(len(newline)):
            f.write(str(int(newline[col]))+' ')
        f.write('\n')
    f.write('\n')
    f.write('Angles\n')
    f.write('\n')
    for nlin in range(len(data.angles.T)):
        newline = data.angles.T[nlin]
        for col in range(len(newline)):
            f.write(str(int(newline[col]))+' ')
        f.write('\n')
    f.write('\n')
    f.write('Dihedrals\n')
    f.write('\n')
    for nlin in range(len(data.dihedrals.T)):
        newline = data.dihedrals.T[nlin]
        for col in range(len(newline)):
            f.write(str(int(newline[col]))+' ')
        f.write('\n')
    f.write('\n')
    f.write('Impropers\n')
    f.write('\n')
    for nlin in range(len(data.impropers.T)):
        newline = data.impropers.T[nlin]
        for col in range(len(newline)):
            f.write(str(int(newline[col]))+' ')
        f.write('\n')
    f.close()
    
