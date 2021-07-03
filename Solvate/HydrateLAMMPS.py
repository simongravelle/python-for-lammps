import numpy as np

def hydrate(solute,solvent,lx,ly,lz,d):
    xcoor = solute.xloxhi
    ycoor = solute.yloyhi
    zcoor = solute.zlozhi

    coordinates = solute.AtomsProperties.coordinates

    Matoms = solute.AtomsProperties.AllAtomsProperties
    nAt = solute.Natoms
    Tat = solute.Tatoms

    Nmolecule = np.max(solute.AtomsProperties.moleculeID)

    Mbonds = solute.Bonds
    nBo = solute.Nbonds
    Tbo = solute.Tbonds

    Mangles = solute.Angles
    Nag = solute.Nangles
    Tag = solute.Tangles

    Mdihedrals = solute.Dihedrals
    Ndi = solute.Ndihedrals
    Tdi = solute.Tdihedrals

    Mimpropers = solute.Impropers
    Nim = solute.Nimpropers
    Tim = solute.Timpropers
    
    addedmol = 0

    molcoordinates = solvent.AtomsProperties.coordinates
    molatoms = solvent.AtomsProperties.AllAtomsProperties
    molatoms[2] += Tat
    Tat += solvent.Tatoms

    molbond = solvent.Bonds
    if type(molbond) != int :
        molbond[1] += Tbo
    Tbo += solvent.Tbonds

    molangle = solvent.Angles
    if type(molangle) != int :
        molangle[1] += Tag
    Tag += solvent.Tangles

    moldihedrals = solvent.Dihedrals
    if type(moldihedrals) != int :
        moldihedrals[1] += Tdi
    Tdi += solvent.Tdihedrals

    molimpropers = solvent.Impropers
    if type(molimpropers) != int :
        molimpropers[1] += Tim
    Tim += solvent.Timpropers
    
    x = -lx / 2
    while x < lx / 2:
        x += d
        y = -ly / 2
        while y < ly / 2:
            y += d
            z = -lz / 2
            while z < lz / 2:
                z += d
                dmin = 1e5
                for lin1 in range(len(coordinates.T)):
                    xi = coordinates[0,lin1]
                    yi = coordinates[1,lin1]
                    zi = coordinates[2,lin1]
                    for lin2 in range(len(molcoordinates.T)):
                        xj = molcoordinates[0,lin2]+x
                        yj = molcoordinates[1,lin2]+y
                        zj = molcoordinates[2,lin2]+z
                        dmin = np.min((dmin, np.sqrt(np.abs(xi-xj)+np.abs(yi-yj)+np.abs(zi-zj))))
                if dmin > d*0.8 :
                    addedmol += 1
                    Nmolecule +=1
                    # add bond
                    if type(molbond) != int :
                        for idx in range(len(molbond.T)) :
                            nBo += 1
                            newline = [nBo, molbond[1,idx], molbond[2,idx]+nAt, molbond[3,idx]+nAt]
                            Mbonds = np.vstack([Mbonds.T,newline]).T
                    # add angle
                    if type(molangle) != int :
                        for idx in range(len(molangle.T)) :
                            Nag += 1
                            newline = [Nag, molangle[1,idx], molangle[2,idx]+nAt, molangle[3,idx]+nAt, molangle[4,idx]+nAt]
                            Mangles = np.vstack([Mangles.T,newline]).T
                    # add dihedrals
                    if type(moldihedrals) != int :
                        for idx in range(len(moldihedrals.T)) :
                            Ndi += 1
                            newline = [Ndi, moldihedrals[1,idx], moldihedrals[2,idx]+nAt, moldihedrals[3,idx]+nAt, moldihedrals[4,idx]+nAt, moldihedrals[5,idx]+nAt]
                            Mdihedrals = np.vstack([Mdihedrals.T,newline]).T
                    # add impropers
                    if type(molimpropers) != int :
                        for idx in range(len(molimpropers.T)) :
                            Nim += 1
                            newline = [Nim, molimpropers[1,idx], molimpropers[2,idx]+nAt, molimpropers[3,idx]+nAt, molimpropers[4,idx]+nAt, molimpropers[5,idx]+nAt]
                            Mimpropers = np.vstack([Mimpropers.T,newline]).T
                    # add atoms
                    if type(molcoordinates) != int :
                        for idx in range(len(molcoordinates.T)) :
                            nAt += 1
                            xj = molcoordinates[0,idx]+x
                            yj = molcoordinates[1,idx]+y
                            zj = molcoordinates[2,idx]+z
                            typ = molatoms[2,idx]
                            q = molatoms[3,idx]
                            newline = [nAt, Nmolecule, typ, q, xj, yj, zj]
                            Matoms = np.vstack([Matoms.T,newline]).T
    class B:
        Natoms = nAt
        Nbonds = nBo
        Nangles = Nag
        Ndihedrals = Ndi
        Nimpropers = Nim
        Tatoms = Tat
        Tbonds = Tbo
        Tangles = Tag
        Tdihedrals = Tdi
        Timpropers = Tim
        atoms  = Matoms
        bonds = Mbonds
        angles = Mangles
        dihedrals = Mdihedrals
        impropers = Mimpropers
        
    LAMMPSData = B()
    return LAMMPSData
