#! /usr/bin/env python
import sys
import os
from numpy import *
from parmed.utils import netcdf as NetCDF

# General class inherited by pdb only for now


class General():
    # converts a string of the type "1,5-8" into a list [1,5,6,7,8]

    def filter2list(self, str):
        residues = [i.split('-') for i in str.split(',')]
        residues = [[int(i) for i in j] for j in residues]
        tmp = []
        for i in residues:
            if len(i) == 1:
                tmp.append(i[0])
            if len(i) == 2:
                tmp.extend(list(range(i[0], i[1] + 1)))
        return tmp


class prmtop:

    def __init__(self, ifile):
        self.file = ifile

#######################
# READS amber prmtop file and extracts specified data.
# Arguments:
#	ifile - name of prmtop file to parse
#	flag - FLAG string identifying data to get (ex: 'FLAG ATOM_NAME')
#	convert - 1=will convert data to numpy array. 0(default) will leave as string list
#######################
    def TopoGet(self, ifile, flag, convert=0):
        f = open(ifile, 'r')
        textfile = f.read()
        f.close()
        n = 0
        for item in textfile.split('%'):
            if n == 1:
                data = item.split()
                data.pop(0)
                break
            if flag in item:
                n = 1
        if convert == 1:
            data = array(list(float(i) for i in data))

        ################################################################
        ## This was another way to do it. Time was about double.       #
        ################################################################
        #~ mass=[]
        #~ while 1:
            #~ x=f.readline()
            #~ if flag in x:
            #~ f.readline()
            #~ while 1:
            #~ y=f.readline()
            #~ if 'FLAG' in y:
            #~ break
            #~ else:
            #~ mass.append(list(float(a) for a in y.split()))
            #~ if not x:
            #~ break
        #~ f.close()
        #~ mass=array(flatten(mass))
        #################################
        return data

    def Get_Masses(self):
        return self.TopoGet(self.file, 'FLAG MASS', 1)

    def Get_Pointers(self):
        return self.TopoGet(self.file, 'FLAG POINTERS')

    def Get_Natoms(self):
        pt = self.Get_Pointers()
        return pt[0]

    def Get_AtomTypes(self):
        return self.TopoGet(self.file, 'FLAG AMBER_ATOM_TYPE')

    def Get_Residues(self):
        return self.TopoGet(self.file, 'FLAG RESIDUE_LABEL')

    def Get_AtomNames(self):
        f = open(self.file, 'r')
        textfile = f.read()
        f.close()
        n = 0
        for item in textfile.split('%'):
            if n == 1:
                x = []
                data = item.split('\n')

                data.pop(0)
                #~ print data
                for line in data:
                    for i in range(20):
                        x.append(line[i * 4:(i + 1) * 4])
                while '' in x:
                    x.remove('')
                break
            if 'FLAG ATOM_NAME' in item:
                n = 1
        x = [i.strip() for i in x]
        return x

    def H_Mask(self):
        atoms = self.Get_AtomNames()
        h_mask = []
        for atomname in atoms:
            if atomname[0] in ['H', 'h']:
                h_mask.append(False)
            else:
                h_mask.append(True)
        return array(h_mask)


#~ A=prmtop('tip4pew.prmtop')
#~ A=prmtop('/home/pjanowsk/York/hairpin/boom/topo.prmtop')
#~ A.masses=A.Get_Masses()
#~ A.natoms=A.Get_Natoms()
#~ print A.masses
#~ print A.file
#~ print A.natoms


class rst7:

    def __init__(self, ifile):
        self.file = ifile
        #~ self.box=self.Get_Box()
        #~ self.natoms=self.Get_Natoms()
        #~ self.box=self.Get_Box()
        #~ self.coords=self.Get_Coords()

    def Get_Natoms(self):
        f = open(self.file, 'r')
        f.readline()
        natoms = f.readline().strip().split()
        natoms = int(natoms[0])
        f.close()
        return natoms

    def Get_Box(self):
        f = open(self.file, 'r')
        x = f.readlines()
        box=x[len(x)-1].split()
        box=array([float(i) for i in box])
        f.close()
        return box

    #~ def Get_Coords(self):
        #~ natoms=self.Get_Natoms()
        #~ with open(self.file) as f:
        #~ x = f.readlines()
        #~ x=x[2:(natoms+1)/2+2]
        #~ y=[]
        #~ for line in x:
        #~ y.extend(line.split())
        #~ y=array([float(coord) for coord in y])
        #~ y=y.reshape(-1,3)
        #~ return y

        ################################################################
        ## This was two other ways to do it. Time was slightly longer.    #
        ################################################################
        #~ def Get_Data(file):
        #~ f=open(file, 'r')
        #~ lines = f.read().splitlines()
        #~ lines.pop(0)
        #~ natoms=int(lines.pop(0)[0:7])
        #~ b=lines.pop()
        #~ b=[b[i:i+12] for i in range(0, len(b), 12)]
        #~ box= array([float(i) for i in b])
        #~ coords=[]
        #~
        #~ for line in range((natoms+1)/2):
        #~ b=[lines[line][i:i+12] for i in range(0, len(lines[line]), 12)]
        #~ coords.append([float(i) for i in b])
        #~ coords=array(coords)
        #~ if natoms%2==1:
        #~ coords=reshape(coords,((natoms+1),3))
        #~ coords=delete(coords,shape(coords)[0]-1,0)
        #~ else:
        #~ coords=reshape(coords,(natoms,3))
        #~ f.close()
        #~ return (natoms,box,coords)
        ###############################################################

    def Get_Coords(self):
        natoms = self.Get_Natoms()
        coords = zeros((natoms, 3))
        f = open(self.file, 'r')
        f.readline()
        f.readline()
        for i in range((natoms + 1) / 2):
            l = float(f.read(12))
            coords[i * 2, 0] = l
            l = float(f.read(12))
            coords[i * 2, 1] = l
            l = float(f.read(12))
            coords[i * 2, 2] = l
            if natoms % 2 == 1 and i == ((natoms + 1) / 2 - 1):
                break
            l = float(f.read(12))
            coords[i * 2 + 1, 0] = l
            l = float(f.read(12))
            coords[i * 2 + 1, 1] = l
            l = float(f.read(12))
            coords[i * 2 + 1, 2] = l
            f.readline()
        f.close()
        return coords
        ################################################################

#~ B=rst7('/home/pjanowsk/York/hairpin/boom/equil30.rst7')
#~ B.natoms=B.Get_Natoms()
#~ B.coords=B.Get_Coords()
#~
#~ print B.natoms
#~ print B.coords


class nc:

    def __init__(self, ifile):
        self.file = ifile

    def Get_Traj(self, skip=1):
        ofile = NetCDF.NetCDFFile(self.file, 'r', mmap=False)
        coords = ofile.variables['coordinates']
        x = coords[::skip, :, :]
        print('Trajectory shape is ')
        print(x.shape)
        ofile.close()
        return x

    def Get_Frames(self):
        ofile = NetCDF.NetCDFFile(self.file, 'r', mmap=False)
        coords = ofile.variables['coordinates']
        frames = coords.shape[0]
        ofile.close()
        return frames

    def Get_CellLens(self, skip=1):
        ofile = NetCDF.NetCDFFile(self.file, 'r', mmap=False)
        celllens = ofile.variables['cell_lengths']
        x = celllens[::skip]
        ofile.close()
        return x

    def Get_CellAngles(self, skip=1):
        ofile = NetCDF.NetCDFFile(self.file, 'r', mmap=False)
        celllens = ofile.variables['cell_angles']
        x = celllens[::skip]
        ofile.close()
        return x

#~ ncfile=nc('md1.nc')
#~ frames=ncfile.Get_Frames()
#~ traj=ncfile.Get_Traj()
#~ print frames
#~ print traj


class pdb(General):

    def __init__(self, ifile):
        self.file = ifile
        self.ATOMlines = []
        with open(self.file) as f:
            for line in f:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    self.ATOMlines.append(line)

    def filter(self, lst, by):
        tmp = []
        if by == 'atom_name':
            for i in self.ATOMlines:
                if i[12:16].strip() in lst:
                    tmp.append(i)
        if by == 'residueNo':
            for i in self.ATOMlines:
                if int(i[22:26].strip()) in lst:
                    tmp.append(i)
        self.ATOMlines = tmp[:]

    def Get_SMTRY(self):
        y = []
        with open(self.file) as f:
            for line in f:
                if "REMARK 290   SMTRY" in line:
                    y.append(line.split())
        if len(y) % 3 != 0:
            print("ERROR: number of SMTRY containing lines is not multiple of 3")
        nsym = len(y) / 3
        nsym = int(nsym)  #  needed for range() statement below(?)
        smtry = zeros((nsym, 3, 3))
        for sym in range(nsym):
            for i in range(3):
                smtry[sym, i, 0] = y[sym * 3 + i][4]
                smtry[sym, i, 1] = y[sym * 3 + i][5]
                smtry[sym, i, 2] = y[sym * 3 + i][6]
        tr = zeros((nsym, 3))
        for sym in range(nsym):
            for i in range(3):
                tr[sym, i] = y[sym * 3 + i][7]
        return smtry, tr

    def Get_Box(self):
        boxinfo = 0
        with open(self.file) as f:
            for line in f:
                if line[0:6] == 'CRYST1':
                    line = line.split()
                    box = line[1:7]
                    box = array([float(i) for i in box])
                    boxinfo = 1
        if boxinfo == 0:
            sys.exit('No CRYST1 record in pdb file. Unable to get box information.')
        return box

    def Get_Coords(self):
        #~ y=[]
        #~ with open(self.file) as f:
            #~ for line in f:
                #~ if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    #~ y.append(line)
        natoms = len(self.ATOMlines)
        coords = zeros((natoms, 3))
        for atom in range(natoms):
            coords[atom, 0] = float(self.ATOMlines[atom][30:38])
            coords[atom, 1] = float(self.ATOMlines[atom][38:46])
            coords[atom, 2] = float(self.ATOMlines[atom][46:54])
        return coords

    def Get_AtomNames(self):
        #~ y=[]
        #~ with open(self.file) as f:
            #~ for line in f:
                #~ if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    #~ y.append(line)
        AtomNames = []
        for atom in self.ATOMlines:
            AtomNames.append(atom[12:16].strip())
        return AtomNames

    def Get_ResidueNumbers(self):
        #~ y=[]
        #~ with open(self.file) as f:
            #~ for line in f:
                #~ if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    #~ y.append(line)
        ResidueNo = []
        for atom in self.ATOMlines:
            ResidueNo.append(int(atom[22:26].strip()))
        return ResidueNo

    def Get_ResidueNames(self):
        ResidueNames = []
        for atom in self.ATOMlines:
            ResidueNames.append(atom[17:20].strip())
        return ResidueNames

    def Get_Bfactors(self):
        #~ y=[]
        #~ with open(self.file) as f:
            #~ for line in f:
                #~ if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    #~ y.append(line)
        Bfactors = []
        for atom in self.ATOMlines:
            Bfactors.append(float(atom[60:66].strip()))
        return Bfactors

#~ pdbfile=pdb('/home/pjanowsk/York/hairpin/2p7e/2P7E.pdb_orig')
#~ coords=pdbfile.Get_Coords()
#~ print coords
#~ smtry,tr =pdbfile.Get_SMTRY()
#~ print smtry
#~ print tr


class evecs:

    def __init__(self, ifile):
        self.file = ifile

    def ReadEvecs(self, verbose=False):
        f = open(self.file, 'r')
        n = f.readline().split()
        n = int(n[n.index('nmodes') + 1])
        natoms = int(f.readline().split()[0])
        w = zeros((n))
        v = zeros((n, natoms))
        n_vec = 0
        n_partial = 0
        width = len(f.readline().strip().split())
        while True:
            line = f.readline().strip().split()
            if "****" in line:
                w[n_vec] = f.readline().strip().split()[1]
                break
        while True:
            line = f.readline().strip().split()
            if not line:
                break
            if "****" in line:
                n_vec += 1
                if (verbose):
                    print(n_vec)
                n_partial = 0
                w[n_vec] = f.readline().strip().split()[1]
            else:
                v[n_vec, n_partial:width + n_partial] = line
                n_partial += width
        f.close()
        self.f = w  # frequencies in cm-1
        self.nm = v.T  # normal modes in columns
        return 0

# x=raf.evecs('evecs7.dat')
# x.ReadEvecs()
# print x.f

#######################################################################
# Calculate Center of Mass of molecule                                #
# Arguments:                                                          #
#     coords: nx3 array of cartesian coordinates of each atom         #
#     masses: 1xn array of mass of each atom. If None, then all       #
#             atoms have equal weight                                 #
#######################################################################


def COM(coords, masses):
    return average(coords, axis=0, weights=masses)
#~ com=COM(B.coords,A.masses)
#~ print com

# Another way to do this is (not having recourse to numpy.average):
#~ m=masses/sum(masses)   #normalize weights so that sum(m)=1
#~ dot(m,coords)    #matrix multiply m x coords


#######################################################################
# Calculate Center of Mass of molecule with vmd                       #
# Arguments:                                                          #
#     pdbfile: name of pdb file                                       #
#                                                                     #
#######################################################################
def COM_pdb(pdbfile):
    os.system('vmd -dispdev text -e COM.tcl -args %s' % pdbfile)
    f = open('COM.txt', 'r')
    com = array([float(i) for i in f.readline().strip().split()])
    f.close()
    return com
#~ com=COM_pdb('1rpg.pdb')
#~ print com


#######################################################################
# Calculate average structure of a netcdf trajectory                  #
# Arguments:                                                          #
#     traj: nframesxnatomsx3 trajectory array (from ncfile)                               #
# Returns:                                                            #
#     avg: natomsx3 array of coordinates of the average structure     #
#######################################################################
def AverageStructure(traj):
    avg = average(traj, axis=0)
    return avg

#~ def AverageStructure(ncfile):
    #~ ofile = NetCDF.NetCDFFile(ncfile, 'a', mmap=False)
    #~ coords=ofile.variables['coordinates']
    #~ avg=average(coords,axis=0)
    #~ return avg
    #~ ofile.close()

#######################################################################
# Calculate b-factors of a netcdf trajectory                          #
# Arguments:                                                          #
#     traj: nframesxnatomsx3 trajectory array (from ncfile)           #
#     avg: natomsx3 array (output of AverageStructure())              #
# Returns:                                                            #
#     y: natomsx1 array of bfactor for each atom                    #
#######################################################################


def BFactors(traj, avg):
    if avg.shape[0] != traj.shape[1]:
        print("Error: number of atoms in average structure and trajectory do not agree")
    if avg.shape[1] != 3:
        print("Error: average structure atoms do not have 3 dimensions")
    frames = traj.shape[0]
    for i in range(3):
        traj[:, :, i] = traj[:, :, i] - avg[:, i]
    traj = square(traj)
    y = sum(sum(traj, axis=0), axis=1)
    #~ y=sum(y,axis=1)
    y = y / frames * 8 / 3 * (pi**2)
    if y.shape[0] != traj.shape[1]:
        print("Error: no. of B-factors does not agree with number of atoms in trajectory")
    return y

#~ avg=AverageStructure(traj)
#~ bfacs=BFactors(traj,avg)
#~ print bfacs
#~ print avg


#######################################################################
# Calculate rmsd between two coordinate arrays                        #
# Arguments:                                                          #
#     str1: natomsx3 array                                            #
#     str2: natomsx3 array 									          #
#                                                                     #
#######################################################################
def RMSD(str1, str2, masses=0):
    assert(str1.shape == str2.shape), "Error: coordinate array sizes do not match"
    x = square(str1 - str2)
    x = sum(x, axis=1)
    if array_equal(masses, 0) == 1:
        x = average(x)
    else:
        x = average(x, weights=masses)
    return sqrt(x)

#######################################################################
#Calculate rmsd between two coordinate arrays using Kabsch optimal    #
#   alignment (from singular values, this rmsd is NOT mass weighed).  #
#	This is only here for curiosity's sake. For practical use I would #
#	stick to using KabschAlign() followed by RMSD().
# Arguments:                                                          #
#     str1: natomsx3 array                                            #
#     str2: natomsx3 array 									          #
#                                                                     #
#######################################################################


def KabschRMSD(str1, str2, masses=0):
    # create masses array (all ones by default)
    if array_equal(masses, 0) == 1:
        masses = ones(str1.shape[0])
    assert(str1.shape[1] == 3)
    assert (str1.shape[0] == masses.shape[0]
            ), "Error: masses and coordinate lenghts do not match"

    # align center of masses to origin
    str1 = str1 - COM(str1, masses)
    str2 = str2 - COM(str2, masses)
    # correlation matrix
    R = dot(transpose(str2), str1)
    print(R)
    # SVD
    v, s, w_tr = linalg.svd(R)
    print(s)
    # check for inversion
    if linalg.det(v) * linalg.det(w_tr) < 0.0:
        s[-1] = -s[-1]
    # calculate transformation matrix
    E0 = sum(sum(str1 * str1)) + sum(sum(str2 * str2))
    rmsd = sqrt((E0 - 2.0 * sum(s)) / float(str1.shape[0]))
    return rmsd

#######################################################################
# Optimally align coordinates 2 to coordinates1 using Kabsch SVD. For #
# usage example see test_rmsd.py                                      #
# Arguments:                                                          #
#     str1: natomsx3 array                                            #
#     str2: natomsx3 array 									          #
#     masses (optional): 1xn array of masses of each atom    		  #
# Returns:                                                            #
#     coordinates of str2 rmsd aligned to str1                        #
#######################################################################


def KabschAlign(str1, str2, masses=0):
    if array_equal(masses, 0) == 1:
        masses = ones(str1.shape[0])
    assert(str1.shape[1] == 3), "Error: coordinates are not 3-dimensional"
    assert(str1.shape == str2.shape), "Error: coordinate array sizes do not match"
    assert (str1.shape[0] == masses.shape[0]
            ), "Error: masses and coordinate lenghts do not match"

    # align center of masses at origin (translation)
    str1_com = COM(str1, masses)
    str1 = str1 - str1_com
    str2 = str2 - COM(str2, masses)
    # correlation matrix
    R = dot(transpose(str2), str1)
    #~ print R
    # SVD
    v, s, w_tr = linalg.svd(R)
    # check for inversion
    if linalg.det(v) * linalg.det(w_tr) < 0.0:
        # print 'hello'
        v[:, -1] = -v[:, -1]
    # calculate rotation matrix
    U = dot(v, w_tr)
    # rotate
    str2_aligned = dot(str2, U)
    # translate back to center of mass of str1
    str2_aligned = str2_aligned + str1_com
    return str2_aligned


################################################################

def printmdcrd(coords_new, boxsizeline, filename, comment='mdcrd written by pawel'):
    f = open(filename, 'w')
    f.write('%s\n' % comment)
    i = 1
    for atom in coords_new:
        for coord in atom:
            f.write('%8.3f' % coord)
            if i == 10:
                f.write('\n')
                i = 0
            i += 1
    f.write('\n' + boxsizeline + '\n')
    f.close()


def printrst7(coords, box, natoms, filename, comment='rst7 written by pawel'):
    f = open(filename, 'w')
    f.write('%s\n' % comment)
    f.write('%5d\n' % natoms)
    i = 1
    for atom in coords:
        for coord in atom:
            f.write('%12.7f' % coord)
            if i == 6:
                f.write('\n')
                i = 0
            i += 1
    for dim in box:
        f.write('%12.7f' % dim)
    f.write('\n')
    f.close()


def printnc(traj, filename):
    return 1


def pdb_line_out(
        record="ATOM  ",
        atm_no=1,
        atm_name=" CA ",
        alt_loc=" ",
        res_name="ALA",
        chain=" ",
        res_no=1,
        insert_code=" ",
        X=1.0,
        Y=1.0,
        Z=1.0,
        occupancy=1.0,
        B=0.0,
        segment_id="    ",
        element=" H",
        charge="  "):
    print(("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s" % (
        record, atm_no, atm_name, alt_loc, res_name, chain, res_no, insert_code,
        X, Y, Z, occupancy, B, segment_id, element, str(charge))))
    return 0


#######################################################################
# Calculate volume from box vectors and angles                        #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          This is usually taken from rst7 file class, Get_Box routine#
#          Angles must be in degrees.								  #
# Returns:                                                            #
#     volume: volume of the box in A^3	                              #
#######################################################################
def Get_volume(box):
    box[3:6] = radians(box[3:6])
    a, b, c, alpha, beta, gamma = box[0], box[1], box[2], box[3], box[4], box[5],
    b = box[1]
    c = box[2]
    alpha = box[3]
    beta = box[4]
    gamma = box[5]
    volume = a * b * c * sqrt(1 - cos(alpha)**2 - cos(beta) **
                              2 - cos(gamma)**2 + 2 * cos(alpha) * cos(beta) * cos(gamma))

    return volume
#~ B=rst7('/home/pjanowsk/York/hairpin/boom/equil30.rst7')
#~ box=B.Get_Box()
#~ volume=Get_volume(box)


#######################################################################
# Calculate deorthogonalization matrix u and orthogonalization matrix #
# 	invu. Invu is used to take fractional coordinates into cartesian  #
#	space and u is used to take cartesian coordinates into fractional   #
#	space.                                                              #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          This is usually taken from rst7 file class, Get_Box routine#
#          Angles must be in degrees.								                  #
# Returns:                                                            #
#     u: 3x3 array                                                    #
#	  invu: 3x3 array	                                                  #
#######################################################################
def CompXfrm(box):
    box = box.astype(float)
    box[3:6] = radians(box[3:6])
    a, b, c, alpha, beta, gamma = box[0], box[1], box[2], box[3], box[4], box[5]
    volume = a * b * c * sqrt(1 - cos(alpha)**2 - cos(beta) **
                              2 - cos(gamma)**2 + 2 * cos(alpha) * cos(beta) * cos(gamma))
    invu = zeros((3, 3))
    invu[0, 0] = a
    invu[0, 1] = b * cos(gamma)
    invu[0, 2] = c * cos(beta)
    invu[1, 1] = b * sin(gamma)
    invu[1, 2] = c * (cos(alpha) - (cos(beta) * cos(gamma))) / sin(gamma)
    invu[2, 2] = volume / (a * b * sin(gamma))
    u = linalg.inv(invu)
    return u, invu
#~ B=rst7('/home/pjanowsk/York/hairpin/boom/equil30.rst7')
#~ box=B.Get_Box()
#~ u,invu=CompXfrm(box)


#######################################################################
# Compute the metric tensor and reciprocal metric tensor from unit    #
#   cell parameters. The metric tensor (or covariant metric tensor)   #
#   converts reciprocal basis vectors to real space (covariant) basis #
#   vectors. It also operates on real space vectors vectors           #
#   (contravariant components) to produce reciprocal space vectors    #
#   (covariant components). gstar is the reciprocal metric tensor (or #
#   contravariant metric tensor). V is the volume of the real space   #
#   unit cell.                                                        #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          Angles must be in degrees.								                  #
# Returns:                                                            #
#     gmetr: 3x3 array. Metric tensor.                                #
#     gstar: 3x3 array. Reciprocal metric tensor.                     #
#     V: float. Volume of real space unit cell.                       #
#######################################################################
def MetricTensor(box):
    box = box.astype(float)
    box[3:6] = radians(box[3:6])
    a, b, c, alpha, beta, gamma = box[0], box[1], box[2], box[3], box[4], box[5]
    gmetr = zeros((3, 3))
    gmetr[0, 0] = a * a
    gmetr[1, 1] = b * b
    gmetr[2, 2] = c * c
    gmetr[0, 1] = gmetr[1, 0] = a * b * cos(gamma)
    gmetr[0, 2] = gmetr[2, 0] = a * c * cos(beta)
    gmetr[1, 2] = gmetr[2, 1] = b * c * cos(alpha)
    gstar = linalg.inv(gmetr)
    V = sqrt(linalg.det(gmetr))
    return gmetr, gstar, V

#######################################################################
# Compute the dot product of two vectors in any basis (could be a non #
#   orthogonal basis.                                                 #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          Angles must be in degrees.								                  #
#     u: 1x3 array.                                                   #
#     v: 1x3 array.                                                   #
# Returns:                                                            #
#     uv: float. Dot product of u and v in basis defined by box.      #
#######################################################################


def DotProduct(box, u, v):
    gmetr, gstar, V = MetricTensor(box)
    uv = 0
    for i in range(3):
        for j in range(3):
            uv += u[i] * v[j] * gmetr[i, j]
    return uv

#######################################################################
# Compute the reciprocal space unit cell parameters.                  #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          Angles must be in degrees.								                  #
# Returns:                                                            #
#     recip_box: 1x6 array. Reciprocal unit cell box. Angles are in   #
#                degrees of real space unit cell.                     #
#######################################################################


def RecipBox(box):
    gmetr, gstar, V = MetricTensor(box)
    recip_box = zeros((6))
    recip_box[0] = sqrt(DotProduct(box, gstar[0:3, 0], gstar[0:3, 0]))
    recip_box[1] = sqrt(DotProduct(box, gstar[0:3, 1], gstar[0:3, 1]))
    recip_box[2] = sqrt(DotProduct(box, gstar[0:3, 2], gstar[0:3, 2]))
    recip_box[3] = arccos(DotProduct(box, gstar[0:3, 1], gstar[0:3, 2]
                                     ) / (recip_box[1] * recip_box[2])) * 180 / pi
    recip_box[4] = arccos(DotProduct(box, gstar[0:3, 0], gstar[0:3, 2]
                                     ) / (recip_box[0] * recip_box[2])) * 180 / pi
    recip_box[5] = arccos(DotProduct(box, gstar[0:3, 0], gstar[0:3, 1]
                                     ) / (recip_box[0] * recip_box[1])) * 180 / pi
    return recip_box

#######################################################################
# Calculate the isotropic B-factor and isotropic U_equiv given  six   #
#   anisotropic ADPs and the unit cell parameters.                    #
# Arguments:                                                          #
#     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
#          Angles must be in degrees.								                  #
#     Uij: 1x6 array of six anisotropic displacement parameters.      #
#          These are the ANISOU numbers divided by 10000.             #
# Returns:                                                            #
#     B: isotropic B-factor                                           #
#     Ueq: isotropic U_equiv                                          #
#######################################################################


def Adp2B(box, Uij):
    recip_box = RecipBox(box)
    Ueq = Uij[0] * box[0] * box[0] * recip_box[0] * recip_box[0] + \
        Uij[1] * box[1] * box[1] * recip_box[1] * recip_box[1] + \
        Uij[2] * box[2] * box[2] * recip_box[2] * recip_box[2] + \
        2 * Uij[5] * box[0] * box[1] * recip_box[0] * recip_box[1] * cos(box[5] * pi / 180) + \
        2 * Uij[4] * box[0] * box[2] * recip_box[0] * recip_box[2] * cos(box[4] * pi / 180) + \
        2 * Uij[3] * box[1] * box[2] * recip_box[1] * \
        recip_box[2] * cos(box[3] * pi / 180)
    Ueq /= 3
    B = Ueq * 8 * pi * pi
    return B, Ueq


#######################################################################
# Calculate density of a crystal                                      #
# Arguments:                                                          #
#     topo: parmtop file with mass information                        #
#     coord: either .pdb file with CRYST record or .rst7 file with 	  #
#            box information on the last line                         #
# Returns:                                                            #
#     density: floating point number                                  #
#######################################################################
def CalcDensity(topo, coord):
    A = prmtop(topo)
    masses = A.Get_Masses()
    print(('Number of atoms in prmtop: %d' % shape(masses)[0]))
    if coord[-3:] == 'pdb':
        B = pdb(coord)
        print(('Number of atoms in pdb file: %d' % (shape(B.Get_Bfactors())[0])))
        box = B.Get_Box()
    elif coord[-4:] == 'rst7':
        B = rst7(coord)
        print(('Number of atoms in restart file: %d' % B.Get_Natoms()))
        box = B.Get_Box()
    else:
        sys.exit('Unrecognized coordinate file extension')
    volume = Get_volume(box)
    TotMass = sum(masses)
    density = TotMass / volume * 10.0 / 6.022
    print('Total mass: %d Daltons' % TotMass)
    print('Volume: %d A^3' % volume)
    print('Density: %.3f g/cm^3' % density)
    return density
