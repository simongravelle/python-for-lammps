import os
import re
import numpy as np

def search(word, sentences):
    return [i for i in sentences if re.search(r'\b%s\b' % word, i)]

def searchvariable(keywords,type_variable,lines):

    '''
    description : search a keyword in the LAMMPS input file and return the value next to it
    author: Simon Gravelle
    webpage : simongravelle.github.io/
    LAMMPS tutorials : https://lammpstutorials.github.io/
    GITHUB : https://github.com/simongravelle
    youtube channel with LAMMPS videos and scripts : youtube.com/channel/UCLmK_9wpyLVpcP7BPgN6BIw
    '''
    
    if type_variable=='single_integer':
        data=search(keywords, lines);
        if len(data)>0:
            variable = int(data[0][0:-1].replace(keywords,''))
        else:
            variable = 0
    if type_variable=='multiple_floats': 
        data=(search(keywords, lines))
        if len(data)>0:
            data=(search(keywords, lines)[0][0:-1].replace(keywords,''));
            variable=re.findall(r'\d+\.\d+', data);
        else:
            variable=0;
            
    return variable
    
def searchmatrix(keywords,lines,nline,ncolumn):

    '''
    description : search a keyword in the LAMMPS input file and return the matrix right below it
    author: Simon Gravelle
    webpage : simongravelle.github.io/
    LAMMPS tutorials : https://lammpstutorials.github.io/
    GITHUB : https://github.com/simongravelle
    youtube channel with LAMMPS videos and scripts : youtube.com/channel/UCLmK_9wpyLVpcP7BPgN6BIw
    '''

    data = search(keywords, lines)
    if len(data)>0:
        data = search(keywords, lines)[0][0:-1]
        for num, line in enumerate(lines, 1):
            if keywords in line:
                if 'BondBond Coeffs' and 'BondAngle Coeffs' and 'AngleAngle Coeffs' not in line:
                    lineId = num;
        if ncolumn == 0:
        	variable = np.zeros((nline,10))
        else :
        	variable = np.zeros((nline,ncolumn))
        i = 0;
        do = 'yes';
        while do == 'yes':
            if lineId+1+i<len(lines):
                temp = lines[lineId+1+i][0:-1].split(' ');
                data = [x for x in temp if x];
                if len(data) > 0:
                    if len(data) > ncolumn: #to remove the three additional columns that sometimes appear
                        variable[i][0:len(data)] = data[0:7]
                        i += 1
                    else:
                        variable[i][0:len(data)] = data 
                        i += 1
                else:
                    do = 'no';
            else:
                do = 'no';
        # remove null lines
        variable = variable[~np.all(variable == 0, axis=1)]
        if ncolumn == 0:
            # remove null columns 
            nullcolumn = np.all(variable != 0, axis=0)
            for idx in range(len(nullcolumn)):
                status = np.all(variable.T[idx]);
                if status == True:
                    col = idx;
                variable = variable.T[0:col+1]
        else : 
            variable = variable.T
    else:
        variable = 0;
    return variable

def readdatalammps(datafile,atomStyle=None):

    '''
    description : Read data file in LAMMPS format, and return all the values inside it
    author: Simon Gravelle
    webpage : simongravelle.github.io/
    LAMMPS tutorials : https://lammpstutorials.github.io/
    GITHUB : https://github.com/simongravelle
    youtube channel with LAMMPS videos and scripts : youtube.com/channel/UCLmK_9wpyLVpcP7BPgN6BIw
    '''

    with open(datafile) as f:
        datalines = f.readlines();

    #  Header section
    keywordsHeader=['atoms', # number of atoms in system
            'bonds', # number of bonds in system
            'angles', # number of angles in system
            'dihedrals', # number of dihedrals in system
            'impropers', # number of impropers in system
            'atom types', # number of atom types in system
            'bond types', # number of bond types in system
            'angle types', # number of angle types in system
            'dihedral types', # number of dihedral types in system
            'improper types', # number of improper types in system
            'extra bond per atom', # leave space for this many new bonds per atom (deprecated, use extra/bond/per/atom keyword)
            'extra angle per atom', # leave space for this many new angles per atom (deprecated, use extra/angle/per/atom keyword)
            'extra dihedral per atom', # leave space for this many new dihedrals per atom (deprecated, use extra/dihedral/per/atom keyword)
            'extra improper per atom', # leave space for this many new impropers per atom (deprecated, use extra/improper/per/atom keyword)
            'extra special per atom', # leave space for this many new special bonds per atom (deprecated, use extra/special/per/atom keyword)
            'ellipsoids', # number of ellipsoids in system
            'lines', # number of line segments in system
            'triangles', # number of triangles in system
            'bodies', # number of bodies in system
            'xlo xhi', # simulation box boundaries in x dimension
            'ylo yhi', # simulation box boundaries in y dimension
            'zlo zhi', # simulation box boundaries in z dimension
            'xy xz yz' # simulation box tilt factors for triclinic system
            ];
    
    keywordsAtomProperty=['Atoms',
            'Velocities',
            'Masses',
            'Ellipsoids',
            'Lines',
            'Triangles',
            'Bodies'
            ];

    keywordsMoleculartopology=['Bonds',
            'Angles',
            'Dihedrals',
            'Impropers',
            ];   
    
    keywordsForcefield=['Pair Coeffs',
            'PairIJ Coeffs',
            'Bond Coeffs',
            'Angle Coeffs',
            'Dihedral Coeffs',
            'Improper Coeffs'
            ];   
        
    class A:
        
        ############
        ## Header ##
        ############
        
        Natoms = searchvariable(keywordsHeader[0],'single_integer',datalines);
        Nbonds = searchvariable(keywordsHeader[1],'single_integer',datalines);
        Nangles = searchvariable(keywordsHeader[2],'single_integer',datalines);
        Ndihedrals = searchvariable(keywordsHeader[3],'single_integer',datalines);
        Nimpropers = searchvariable(keywordsHeader[4],'single_integer',datalines);
        Tatoms = searchvariable(keywordsHeader[5],'single_integer',datalines);
        Tbonds = searchvariable(keywordsHeader[6],'single_integer',datalines);
        Tangles = searchvariable(keywordsHeader[7],'single_integer',datalines);
        Tdihedrals = searchvariable(keywordsHeader[8],'single_integer',datalines);
        Timpropers = searchvariable(keywordsHeader[9],'single_integer',datalines);
        global nline 
        nline = np.max([Natoms,Nbonds,Nangles,Ndihedrals,Nimpropers,Tatoms,Tbonds,Tangles,Tdihedrals,Timpropers])
        extraBond = searchvariable(keywordsHeader[10],'single_integer',datalines);
        extraAngle = searchvariable(keywordsHeader[11],'single_integer',datalines);
        extraDihedral = searchvariable(keywordsHeader[12],'single_integer',datalines);
        extraImproper = searchvariable(keywordsHeader[13],'single_integer',datalines);
        extraSpecial = searchvariable(keywordsHeader[14],'single_integer',datalines);
        ellipsoids = searchvariable(keywordsHeader[15],'single_integer',datalines);
        lines = searchvariable(keywordsHeader[16],'single_integer',datalines);
        triangles = searchvariable(keywordsHeader[17],'single_integer',datalines);
        bodies = searchvariable(keywordsHeader[18],'single_integer',datalines);
        xloxhi = searchvariable(keywordsHeader[19],'multiple_floats',datalines);
        yloyhi = searchvariable(keywordsHeader[20],'multiple_floats',datalines);
        zlozhi = searchvariable(keywordsHeader[21],'multiple_floats',datalines);
        xyxzyz = searchvariable(keywordsHeader[22],'multiple_floats',datalines);
        
        ############################
        ## Atom property sections ##
        ############################
        
        class B:
            AllAtomsProperties = searchmatrix(keywordsAtomProperty[0],datalines,nline,7)
            if atomStyle == 'full' or atomStyle == None:
                if atomStyle == None:
                    print('Atom style full has been assumed')
                atomID = AllAtomsProperties[0]
                moleculeID = AllAtomsProperties[1]
                atomtype = AllAtomsProperties[2]
                charge = AllAtomsProperties[3]
                coordinates = AllAtomsProperties[4:7]
            elif atomStyle == 'atomic':
                atomID = AllAtomsProperties[0]
                atomtype = AllAtomsProperties[1]
                coordinates = AllAtomsProperties[2:5] 
            elif atomStyle == 'molecular':
                atomID = AllAtomsProperties[0]
                moleculeID = AllAtomsProperties[1]
                atomtype = AllAtomsProperties[2]
                coordinates = AllAtomsProperties[3:6]
        AtomsProperties = B()
        Velocities = searchmatrix(keywordsAtomProperty[1],datalines,nline,4);
        Masses = searchmatrix(keywordsAtomProperty[2],datalines,nline,2);
        # Still to be implemented : Ellipsoids, Lines, Triangles, Bodies
        
        #################################
        ## Molecular topology sections ##
        #################################
        Bonds = searchmatrix(keywordsMoleculartopology[0],datalines,nline,4);
        Angles = searchmatrix(keywordsMoleculartopology[1],datalines,nline,5);
        Dihedrals = searchmatrix(keywordsMoleculartopology[2],datalines,nline,6);
        Impropers = searchmatrix(keywordsMoleculartopology[3],datalines,nline,6);

        ##########################
        ## Force field sections ##
        ##########################   
        PairCoeffs = searchmatrix(keywordsForcefield[0],datalines,nline,0);
        PairIJCoeffs = searchmatrix(keywordsForcefield[1],datalines,nline,0);
        AngleCoeffs = searchmatrix(keywordsForcefield[2],datalines,nline,0);
        DihedralCoeffs = searchmatrix(keywordsForcefield[3],datalines,nline,0);
        ImproperCoeffs = searchmatrix(keywordsForcefield[4],datalines,nline,0);
        
        # Still to be implemented :  class 2 force field sections
        # class 2 force field sections
        #BondBond Coeffs, BondAngle Coeffs, MiddleBondTorsion Coeffs, EndBondTorsion Coeffs, AngleTorsion Coeffs, AngleAngleTorsion Coeffs, BondBond13 Coeffs, AngleAngle Coeffs

    LAMMPSData = A()

    return LAMMPSData
