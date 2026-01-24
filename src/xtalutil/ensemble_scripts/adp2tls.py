#! /usr/bin/env python2

#=====================================================================#
# Calculate TLS parameters given a PDB file with Uij or B-factors or  #
#   a covariance matrix from cpptraj. TLS parameters are derived from #
#   a least squares fit to mean square atomic displacements via eq(5) #
#   in Winn et al, "Use of TLS parameters to model anisotropic        #
#   displacements in macromolecular refinement", Acta Cryst D (2001), #
#   D57, 122 which can be expanded to eq(1) in Painter & Merrit, Acta #
#   Cryst. D. (2006), D62, 439. If a covariance matrix file is        #
#   provided, an (6n)*2 set of equations is solved. If only the PDB   #
#   file is given, a 6N set of equations is constructed from          #
#   the ANISOU anisotropic displacement parameters. If ANISOU records #
#   are not present in the file, equivalent ADPs calculated from the  #
#   isotropic B-factor are used.                                      #
#                                                                     #
# This script draws heavily on code found in the Python Macromolecular#
# Library (mmLib), J. Painter and E.A. Merritt, "mmLib Python toolkit #
# for manipulating annotated structural models of biological          #
# macromolecules." J. Appl. Cryst. (2004) 37, 174-178.                #
# In particular, the code for calculating the center of reaction      #
# (center_of_reaction_shift()) is taken from mmLib and only slightly  #
# adapted to interface with this code. We are deeply indebted to the  #
# excellent mmLib library which can be found at                       #
#                         http://pymmlib.sourceforge.net/             #
# and is licensed under The Artistic License 2.0.                     #
#                                                                     #
# Arguments:                                                          #
#     -i input pdb file (if ANISOU present will use that, otherwise   #
#        uses B-factors)                                              #
#     -o prefix of output files. Prints $prefix.pdb and $prefix.dat.  #
#     -cov <filename> (optional) Covariance matrix file.              #
#     -nocor (optional) Skip translation to center of reaction        #
#     -biso (optional) fit to TLS + Biso model, as described in       #
#        Painter & Merritt, Acta Cryst. D62: 439, 2006.               # 
#                                                                     #
# Returns:                                                            #
#     Prints the least squares residual between the original set of   #
#     anisotropic ADP and the ADP calculated from the best-fit TLS    #
#     model. Creates a pdb file ($prefix.pdb) with the TLS-fit        #
#     anisotropic ADPs and a REFMAC style TLS data file               #
#     ($prefix.refmac.txt).                                           #
#                                                                     #
# Example usage:                                                      #
#     adp2tls.py -i 4lzt.pdb -o tls1                                  #
#                                                                     #
# Pawel Janowski, David Case Group, Rutgers U., Jan. 2015             #
#######################################################################

import sys
import os
import numpy as np
from numpy import linalg
import argparse
# from parmed.structure import Structure
from parmed import read_PDB, write_PDB


def MetricTensor(box):
    '''
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
    #          Angles must be in degrees.                                 #
    # Returns:                                                            #
    #     gmetr: 3x3 array. Metric tensor.                                #
    #     gstar: 3x3 array. Reciprocal metric tensor.                     #
    #     V: float. Volume of real space unit cell.                       #
    #######################################################################
    '''
    box = box.astype(float)
    box[3:6] = np.radians(box[3:6])
    a, b, c, alpha, beta, gamma = box[0], box[1], box[2], box[3], box[4], box[5]
    gmetr = np.zeros((3, 3))
    gmetr[0, 0] = a * a
    gmetr[1, 1] = b * b
    gmetr[2, 2] = c * c
    gmetr[0, 1] = gmetr[1, 0] = a * b * np.cos(gamma)
    gmetr[0, 2] = gmetr[2, 0] = a * c * np.cos(beta)
    gmetr[1, 2] = gmetr[2, 1] = b * c * np.cos(alpha)
    gstar = np.linalg.inv(gmetr)
    V = np.sqrt(np.linalg.det(gmetr))
    return gmetr, gstar, V


def DotProduct(box, u, v):
    '''
    #######################################################################
    # Compute the dot product of two vectors in any basis (could be a non #
    #   orthogonal basis.                                                 #
    # Arguments:                                                          #
    #     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
    #          Angles must be in degrees.                                 #
    #     u: 1x3 array.                                                   #
    #     v: 1x3 array.                                                   #
    # Returns:                                                            #
    #     uv: float. Dot product of u and v in basis defined by box.      #
    #######################################################################
    '''
    gmetr, gstar, V = MetricTensor(box)
    uv = 0
    for i in range(3):
        for j in range(3):
            uv += u[i] * v[j] * gmetr[i, j]
    return uv


def RecipBox(box):
    '''
    #######################################################################
    # Compute the reciprocal space unit cell parameters.                  #
    # Arguments:                                                          #
    #     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
    #          Angles must be in degrees.                                 #
    # Returns:                                                            #
    #     recip_box: 1x6 array. Reciprocal unit cell box. Angles are in   #
    #                degrees of real space unit cell.                     #
    #######################################################################
    '''
    gmetr, gstar, V = MetricTensor(box)
    recip_box = np.zeros((6))
    recip_box[0] = np.sqrt(DotProduct(box, gstar[0:3, 0], gstar[0:3, 0]))
    recip_box[1] = np.sqrt(DotProduct(box, gstar[0:3, 1], gstar[0:3, 1]))
    recip_box[2] = np.sqrt(DotProduct(box, gstar[0:3, 2], gstar[0:3, 2]))
    recip_box[3] = np.arccos(DotProduct(box, gstar[0:3, 1], gstar[0:3, 2]) /
                             (recip_box[1] * recip_box[2])) * 180 / np.pi
    recip_box[4] = np.arccos(DotProduct(box, gstar[0:3, 0], gstar[0:3, 2]) /
                             (recip_box[0] * recip_box[2])) * 180 / np.pi
    recip_box[5] = np.arccos(DotProduct(box, gstar[0:3, 0], gstar[0:3, 1]) /
                             (recip_box[0] * recip_box[1])) * 180 / np.pi
    return recip_box


def Adp2B(box, Uij, method=2):
    '''
    #######################################################################
    # Calculate the isotropic B-factor and isotropic U_equiv given  six   #
    #   anisotropic ADPs and the unit cell parameters.                    #
    # Arguments:                                                          #
    #     box: 1x6 array of box vectors [a,b,c,alpha,beta,gamma]          #
    #          Angles must be in degrees.                                 #
    #     Uij: 1x6 array of six anisotropic displacement parameters.      #
    #          These are the ANISOU numbers divided by 10000.             #
    #     method: 1 - used for cif-style Uij (dimensionless ADP in unit   #
    #                 cell basis)                                         #
    #             2 - used for ANISOU style Uij(cartesian basis           #
    #                 displacement in A**2) (default)                     #
    # Returns:                                                            #
    #     B: isotropic B-factor                                           #
    #     Ueq: isotropic U_equiv                                          #
    #######################################################################
    '''
    if method == 2:
        Ueq = (Uij[0] + Uij[1] + Uij[2]) / 3.
    else:
        recip_box = RecipBox(box)
        Ueq = Uij[0] * box[0] * box[0] * recip_box[0] * recip_box[0] + \
            Uij[1] * box[1] * box[1] * recip_box[1] * recip_box[1] + \
            Uij[2] * box[2] * box[2] * recip_box[2] * recip_box[2] + \
            2 * Uij[5] * box[0] * box[1] * recip_box[0] * recip_box[1] * np.cos(box[5] * np.pi / 180) + \
            2 * Uij[4] * box[0] * box[2] * recip_box[0] * recip_box[2] * np.cos(box[4] * np.pi / 180) + \
            2 * Uij[3] * box[1] * box[2] * recip_box[1] * \
            recip_box[2] * np.cos(box[3] * np.pi / 180)
        Ueq /= 3
    B = Ueq * 8 * np.pi * np.pi
    return B, Ueq


def B2Adp(Bfactor):
    """
    Convert isotropic ADP to anisotropic ADP equivalent.
    """
    Uij = np.zeros((6))
    U = Bfactor * (3. / 8) / (np.pi * np.pi) / 3.
    Uij[0:3] = np.array([U, U, U])
    return Uij


def center_of_mass(pdb, masses=None):
    """
    #######################################################################
    # Calculate center of atoms in the pdb structure. If masses is None   #
    # all atoms have equal weight. Optionally provide 1xN array of masses #
    # to get the mass-weighed center of mass                              #
    #######################################################################
    """
    natoms = len(pdb.atoms)
    xyz = pdb.coordinates
    com = np.average(xyz, axis=0, weights=masses)
    return com


def set_A(pdb, com, w=1):
    '''
    A is a 6Nx20 matrix. For each atom, sets six rows of matrix A with
    the TLS coefficients for an atom located at position x,y,z with
    least-squares weight w. This is according to eq 5 in Winn et al.,
    "Use of TLS parameters...", Acta Cryst D, 2000. eq 1 in Painter and
    Merritt, "Optimal description of a protein structure...", Acta Cryst
    D, 2006

    If biso is set, add N extra columns to handle the Biso variables.
    '''
    # create zeros matrix
    natoms = len(pdb.atoms)
    xyz = pdb.coordinates
    xyz = xyz - com
    if args.biso:
        A = np.zeros((natoms * 6, 20 + natoms))
    else:
        A = np.zeros((natoms * 6, 20))

    # use label indexing to avoid confusion!
    T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
        S1133, S2211, S12, S13, S23, S21, S31, S32 = \
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)

    for i in range(natoms):
        # coordinates
        x = xyz[i, 0]
        y = xyz[i, 1]
        z = xyz[i, 2]
        xx = x * x
        yy = y * y
        zz = z * z
        xy = x * y
        xz = x * z
        yz = y * z

        # indices of the components of U
        U11 = i * 6
        U22 = i * 6 + 1
        U33 = i * 6 + 2
        U12 = i * 6 + 3
        U13 = i * 6 + 4
        U23 = i * 6 + 5

        # populate A
        A[U11, T11] = w * 1.0
        A[U11, L22] = w * zz
        A[U11, L33] = w * yy
        A[U11, L23] = w * -2.0 * yz
        A[U11, S31] = w * -2.0 * y
        A[U11, S21] = w * 2.0 * z

        A[U22, T22] = w * 1.0
        A[U22, L11] = w * zz
        A[U22, L33] = w * xx
        A[U22, L13] = w * -2.0 * xz
        A[U22, S12] = w * -2.0 * z
        A[U22, S32] = w * 2.0 * x

        A[U33, T33] = w * 1.0
        A[U33, L11] = w * yy
        A[U33, L22] = w * xx
        A[U33, L12] = w * -2.0 * xy
        A[U33, S23] = w * -2.0 * x
        A[U33, S13] = w * 2.0 * y

        A[U12, T12] = w * 1.0
        A[U12, L33] = w * -xy
        A[U12, L23] = w * xz
        A[U12, L13] = w * yz
        A[U12, L12] = w * -zz
        A[U12, S2211] = w * z
        A[U12, S31] = w * x
        A[U12, S32] = w * -y

        A[U13, T13] = w * 1.0
        A[U13, L22] = w * -xz
        A[U13, L23] = w * xy
        A[U13, L13] = w * -yy
        A[U13, L12] = w * yz
        A[U13, S1133] = w * y
        A[U13, S23] = w * z
        A[U13, S21] = w * -x

        A[U23, T23] = w * 1.0
        A[U23, L11] = w * -yz
        A[U23, L23] = w * -xx
        A[U23, L13] = w * xy
        A[U23, L12] = w * xz
        A[U23, S2211] = w * -x
        A[U23, S1133] = w * -x
        A[U23, S12] = w * y
        A[U23, S13] = w * -z

        if args.biso:
            A[U11, i+20] = w
            A[U22, i+20] = w
            A[U33, i+20] = w

    return A


def set_A_cov(pdb, com, w=1):
    """
    A is a (9N(N+1)/2 - 3*N) x21 matrix. Each row corresponds to one <U_iU_j>.
    The covariance matrix consists of NxN 3x3 matrices. Each 3x3 matrix has
    9 unique <U_iU_j>. The N 3x3 matrices on the diagonal have only 6
    unique <U_iU_i> so need to subtract 3*N rows. So in total there are
    [N(N+1)/2] unique 3x3 matrices which gives 9*[N(N+1)/2] - 3*N rows. So
    A is 9*[N(N+1)/2] - 3*N row by 21 columns (21 TLS parameters). For
    each atom pair i and j 9 equations but if i=j 6 equations.

    If biso is set, add N extra columns to handle the Biso variables.
    """
    # create zeros matrix
    N = len(pdb.atoms)
    xyz = pdb.coordinates
    xyz = xyz - com
    if args.biso:
        A = np.zeros((9 * N * (N + 1) / 2 - 3 * N, 21 + N))
    else:
        A = np.zeros((9 * N * (N + 1) / 2 - 3 * N, 21))

    # use label indexing to avoid confusion!
    T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
        S11, S22, S33, S12, S13, S23, S21, S31, S32 = \
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)

    row = 0
    for i in range(N):
        for j in range(i, N):
            # coordinates
            x1, y1, z1 = xyz[i, 0], xyz[i, 1], xyz[i, 2]
            x2, y2, z2 = xyz[j, 0], xyz[j, 1], xyz[j, 2]

            # indices of the components of U
            U11 = row
            U22 = row + 1
            U33 = row + 2
            U12 = row + 3
            U13 = row + 4
            U23 = row + 5
            if i == j:
                row += 6
            else:
                U21 = row + 6
                U31 = row + 7
                U32 = row + 8
                row += 9

            # populate A
            A[U11, T11] = w
            A[U11, L22] = w * z1 * z2
            A[U11, L33] = w * y1 * y2
            A[U11, L23] = w * -1.0 * (y1 * z2 + z1 * y2)
            A[U11, S31] = w * -1.0 * (y1 + y2)
            A[U11, S21] = w * (z1 + z2)

            A[U22, T22] = w
            A[U22, L11] = w * z1 * z2
            A[U22, L33] = w * x1 * x2
            A[U22, L13] = w * -1.0 * (x1 * z2 + z1 * x2)
            A[U22, S12] = w * -1.0 * (z1 + z2)
            A[U22, S32] = w * (x1 + x2)

            A[U33, T33] = w
            A[U33, L11] = w * y1 * y2
            A[U33, L22] = w * x1 * x2
            A[U33, L12] = w * -1.0 * (x1 * y2 + y1 * x2)
            A[U33, S23] = w * -1.0 * (x1 + x2)
            A[U33, S13] = w * (y1 + y2)

            A[U12, T12] = w
            A[U12, L33] = w * -1.0 * y1 * x2
            A[U12, L23] = w * z1 * x2
            A[U12, L13] = w * y1 * z2
            A[U12, L12] = w * -1.0 * z1 * z2
            A[U12, S11] = w * -1.0 * z2
            A[U12, S22] = w * z1
            A[U12, S31] = w * x2
            A[U12, S32] = w * -1.0 * y1

            A[U13, T13] = w
            A[U13, L22] = w * -1.0 * z1 * x2
            A[U13, L23] = w * y1 * x2
            A[U13, L13] = w * -1.0 * y1 * y2
            A[U13, L12] = w * z1 * y2
            A[U13, S11] = w * y2
            A[U13, S33] = w * -1.0 * y1
            A[U13, S23] = w * z1
            A[U13, S21] = w * -1.0 * x2

            A[U23, T23] = w
            A[U23, L11] = w * -1.0 * z1 * y2
            A[U23, L23] = w * -1.0 * x1 * x2
            A[U23, L13] = w * x1 * y2
            A[U23, L12] = w * z1 * x2
            A[U23, S22] = w * -1.0 * x2
            A[U23, S33] = w * x1
            A[U23, S12] = w * y2
            A[U23, S13] = w * -1.0 * z1

            if i != j:
                A[U21, T12] = w
                A[U21, L33] = w * -1.0 * x1 * y2
                A[U21, L23] = w * x1 * z2
                A[U21, L13] = w * z1 * y2
                A[U21, L12] = w * -1.0 * z1 * z2
                A[U21, S11] = w * -1.0 * z1
                A[U21, S22] = w * z2
                A[U21, S31] = w * x1
                A[U21, S32] = w * -1.0 * y2

                A[U31, T13] = w
                A[U31, L22] = w * -1.0 * x1 * z2
                A[U31, L23] = w * x1 * y2
                A[U31, L13] = w * -1.0 * y1 * y2
                A[U31, L12] = w * y1 * z2
                A[U31, S11] = w * y1
                A[U31, S33] = w * -1.0 * y2
                A[U31, S23] = w * z2
                A[U31, S21] = w * -1.0 * x1

                A[U32, T23] = w
                A[U32, L11] = w * -1.0 * y1 * z2
                A[U32, L23] = w * -1.0 * x1 * x2
                A[U32, L13] = w * y1 * x2
                A[U32, L12] = w * x1 * z2
                A[U32, S22] = w * -1.0 * x1
                A[U32, S33] = w * x2
                A[U32, S12] = w * y1
                A[U32, S13] = w * -1.0 * z2
            else:
                if args.biso:
                    A[U11, i+21] = w
                    A[U22, i+21] = w
                    A[U33, i+21] = w

    return A


def set_B(pdb, w=1):
    """
    B has 6*N rows and a single column. For each atom, sets the six
    rows of vector b with the experimental/target anisotropic ADP values
    U with weight w.
    """
    natoms = len(pdb.atoms)
    B = np.zeros((natoms * 6))
    for i, atom in enumerate(pdb.atoms):
        if atom.anisou is None:
            atom.anisou = B2Adp(atom.bfactor)
        B[i * 6:i * 6 + 6] = atom.anisou
    return B


def set_B_cov(pdb, covariance_file, w=1):
    """
    B is 9*[N(N+1)/2] - 3*N rows by 1 column. That's because there are 9
    unique <UiUj> for each atom pair but only 6 when i=j.
    """
    N = len(pdb.atoms)
    B = np.zeros((9 * N * (N + 1) / 2 - 3 * N))
    uu = np.genfromtxt(covariance_file)
    assert uu.shape == (3 * N, 3 * N), "Covariance matrix shape does not match " \
        "number of atoms. Should be 3Nx3N."
    row = 0
    for i in range(N):
        for j in range(i, N):
            B[row] = uu[i * 3, j * 3]
            B[row + 1] = uu[i * 3 + 1, j * 3 + 1]
            B[row + 2] = uu[i * 3 + 2, j * 3 + 2]
            B[row + 3] = uu[i * 3, j * 3 + 1]
            B[row + 4] = uu[i * 3, j * 3 + 2]
            B[row + 5] = uu[i * 3 + 1, j * 3 + 2]
            if i == j:
                row += 6
            else:
                B[row + 6] = uu[i * 3 + 1, j * 3]
                B[row + 7] = uu[i * 3 + 2, j * 3]
                B[row + 8] = uu[i * 3 + 2, j * 3 + 1]
                row += 9
    return B, uu


def solve_SVD(A, B):
    """
    Solve Ax=B for x.
    """
    # solve by SVD
    U, W, Vt = np.linalg.svd(A, full_matrices=0)

    V = np.transpose(Vt)
    Ut = np.transpose(U)

    # analyze singular values and generate smallness cutoff
    cutoff = max(W) * 1E-10

    # make W
    dim_W = len(W)
    Wi = np.zeros((dim_W, dim_W), float)

    for i in range(dim_W):
        if W[i] > cutoff:
            Wi[i, i] = 1.0 / W[i]

    # solve for x: really inefficient here!
    UtB = np.dot(Ut, B)
    WUtB = np.dot(Wi, UtB)
    x = np.dot(V, WUtB)
    return x


def set_TLS(x):
    """
    Convert the 20x1 or 21x1 array of TLS parameters obtained from solving
    the least squares problem and convert to 3x3 T,L,S arrays.
    """

    if args.covariance is None:
        # use label indexing to avoid confusion!
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
            S1133, S2211, S12, S13, S23, S21, S31, S32 = \
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)

        T = np.array([[x[T11], x[T12], x[T13]],
                      [x[T12], x[T22], x[T23]],
                      [x[T13], x[T23], x[T33]]])

        L = np.array([[x[L11], x[L12], x[L13]],
                      [x[L12], x[L22], x[L23]],
                      [x[L13], x[L23], x[L33]]])

        # S11+S22+S33 = 0
        S22 = 2.0 * x[S2211] / 3.0 + x[S1133] / 3.0
        S11 = S22 - x[S2211]
        S33 = S11 - x[S1133]
        S = np.array([[S11, x[S12], x[S13]],
                      [x[S21], S22, x[S23]],
                      [x[S31], x[S32], S33]])

    else:
        # use label indexing to avoid confusion!
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
            S11, S22, S33, S12, S13, S23, S21, S31, S32 = \
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)

        T = np.array([[x[T11], x[T12], x[T13]],
                      [x[T12], x[T22], x[T23]],
                      [x[T13], x[T23], x[T33]]])

        L = np.array([[x[L11], x[L12], x[L13]],
                      [x[L12], x[L22], x[L23]],
                      [x[L13], x[L23], x[L33]]])

        S = np.array([[x[S11], x[S12], x[S13]],
                      [x[S21], x[S22], x[S23]],
                      [x[S31], x[S32], x[S33]]])

    return T, L, S


def write_TLS(T, L, S, com, pdb, prefix):
    '''
    Save T, L, S to file in REFMAC format
    '''
    L = L * (180 / np.pi)**2
    S = S * (180 / np.pi)
    start_res = pdb.atoms[0].residue.number
    start_chain = pdb.atoms[0].residue.chain
    end_res = pdb.atoms[-1].residue.number
    end_chain = pdb.atoms[-1].residue.chain
    f = open("%s.refmac.txt" % prefix, 'w')
    f.write('REFMAC\n\n')
    f.write('TLS\n')
    f.write('RANGE  \'%s%5d.\' \'%s%5d.\' ALL\n'
            % (start_chain, start_res, end_chain, end_res))
    f.write('ORIGIN   %8.4f %8.4f %8.4f\n' % (com[0], com[1], com[2]))
    f.write('T   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' %
            (T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[0, 2], T[1, 2]))
    f.write('L   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' %
            (L[0, 0], L[1, 1], L[2, 2], L[0, 1], L[0, 2], L[1, 2]))
    # <S22 - S11> <S11 - S33> <S12> <S13> <S23> <S21> <S31> <S32>
    f.write('S   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n'
            % (S[1, 1] - S[0, 0], S[0, 0] - S[2, 2], S[0, 1], S[0, 2], S[1, 2], S[1, 0], S[2, 0], S[2, 1]))
    f.close()
    return 0

def write_TLS_remark3(T, L, S, com, pdb, prefix):
    '''
    Save T, L, S to file in PDB-style REMARK 3 format
    '''
    L = L * (180 / np.pi)**2
    S = S * (180 / np.pi)
    start_res = pdb.atoms[0].residue.number
    start_chain = pdb.atoms[0].residue.chain
    if start_chain == "":
        start_chain = "A"
    end_res = pdb.atoms[-1].residue.number
    end_chain = pdb.atoms[-1].residue.chain
    if end_chain == "":
        end_chain = "A"
    f = open("%s.remark3.dat" % prefix, 'w')
    f.write('REMARK   3\n')
    f.write('REMARK   3  TLS DETAILS\n')
    f.write('REMARK   3   NUMBER OF TLS GROUPS  :    1\n')
    f.write('REMARK   3\n')
    f.write('REMARK   3   TLS GROUP :     1\n')
    f.write('REMARK   3    NUMBER OF COMPONENTS GROUP :    1\n')
    f.write('REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI\n')
    f.write('REMARK   3    RESIDUE RANGE :   %1s %5d        %1s %5d\n' %
                              (start_chain, start_res, end_chain, end_res))
    f.write('REMARK   3    ORIGIN FOR THE GROUP (A):%9.4f%9.4f%9.4f\n' %
                              (com[0], com[1], com[2]))
    f.write('REMARK   3    T TENSOR\n')
    f.write('REMARK   3      T11: %8.4f T22: %8.4f\n' % (T[0, 0], T[1, 1]))
    f.write('REMARK   3      T33: %8.4f T12: %8.4f\n' % (T[2, 2], T[0, 1]))
    f.write('REMARK   3      T13: %8.4f T23: %8.4f\n' % (T[0, 2], T[1, 2]))
    f.write('REMARK   3    L TENSOR\n')
    f.write('REMARK   3      L11: %8.4f L22: %8.4f\n' % (L[0, 0], L[1, 1]))
    f.write('REMARK   3      L33: %8.4f L12: %8.4f\n' % (L[2, 2], L[0, 1]))
    f.write('REMARK   3      L13: %8.4f L23: %8.4f\n' % (L[0, 2], L[1, 2]))
    f.write('REMARK   3    S TENSOR\n')
    f.write('REMARK   3      S11: %8.4f S12: %8.4f S13: %8.4f\n' % 
                                 (S[0, 0], S[0, 1], S[0, 2]))
    f.write('REMARK   3      S21: %8.4f S22: %8.4f S23: %8.4f\n' % 
                                 (S[1, 0], S[1, 1], S[1, 2]))
    f.write('REMARK   3      S31: %8.4f S32: %8.4f S33: %8.4f\n' % 
                                 (S[2, 0], S[2, 1], S[2, 2]))
    f.write('REMARK   3\n')
    f.close()
    return 0

def calc_ADP_from_TLS(x, A, B, natoms):
    """
    Back-calculate the ADP predicted by the best-fit TLS parameters and get
    the residual to the original data. UTLS is a 6Nx1 array of the calculated
    Uij, 6 for each atom.
    """
    UTLS = np.dot(A, x)
    SSE = np.sum((UTLS - B)**2)
    MSE = SSE / len(B)
    # pearson correlation coefficinet
    r2 = (np.dot(B,UTLS))**2/(np.sum(B**2) * np.sum(UTLS**2))

    for i in range(len(B)):
        print " %12.5f  %12.5f" % (B[i], UTLS[i])

    if args.biso:
        # follow the procedure in Section 3.1.2 of Painter & Merritt,
        # 2006, to make all the Biso values positive
        #
        #  Note: this is just a local calculation here, for printing to
        #  the ouptut file: it doesn't affect the TLS matrices or ADP's
        #  that get written to the output pdb file.
        #
        #  first, need to get the "pure" TLS motion by setting all
        #   columns of a beyond 20 or 21 to zero:
        if args.covariance:
            startcol = 21
        else:
            startcol = 20
        for row in range(len(B)):
            for column in range(natoms):
                A[row,column+startcol] = 0
        UTLSp = np.dot(A,x)

        print "#   i     Uobs     Utls       Uiso"
        Uisomin = 9999.
        for i in range(natoms):
            Utls = (UTLSp[6*i] + UTLSp[6*i+1] + UTLSp[6*i+2])/3.
            Uobs = (B[6*i] + B[6*i+1] + B[6*i+2])/3.
            Uiso = Uobs - Utls
            print "%5d  %8.3f %8.3f %8.3f" % ( i, Uobs, Utls, Uiso )
            if Uiso < Uisomin:
                Uisomin = Uiso
        print "Uisomin is %12.5f" % Uisomin

    return UTLS, SSE, MSE, r2


def calc_ADP_from_TLS_cov(pdb, com, uu, xv, w=1.0):
    """
    Back-calculate the ADP predicted by the best-fit TLS parameters and get
    the residual to the original data. UTLS is a 6Nx1 array of the calculated
    Uij, 6 for each atom. This is more tricky in the full covariance case. We
    are only looking for contribution from the variance equations (the 6
    ADP for each atom) to make this more comparable to the typical TLS. In
    this case the A matrix has 21 columns and there are 21 TLS parameters
    in x (because full covariance removes degeneracy around S11,S22,S33. So
    we need to rebuild A, B from scratch.
    """

    # build A
    # create zeros matrix
    natoms = len(pdb.atoms)
    xyz = pdb.coordinates
    xyz = xyz - com
    if args.biso:
        A = np.zeros((natoms * 6, 21+natoms))
    else:
        A = np.zeros((natoms * 6, 21))

    # use label indexing to avoid confusion!
    T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
        S11, S22, S33, S12, S13, S23, S21, S31, S32 = \
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)

    for i in range(natoms):
        # coordinates
        x = xyz[i, 0]
        y = xyz[i, 1]
        z = xyz[i, 2]
        xx = x * x
        yy = y * y
        zz = z * z
        xy = x * y
        xz = x * z
        yz = y * z

        # indices of the components of U
        U11 = i * 6
        U22 = i * 6 + 1
        U33 = i * 6 + 2
        U12 = i * 6 + 3
        U13 = i * 6 + 4
        U23 = i * 6 + 5

        # populate A
        A[U11, T11] = w * 1.0
        A[U11, L22] = w * zz
        A[U11, L33] = w * yy
        A[U11, L23] = w * -2.0 * yz
        A[U11, S31] = w * -2.0 * y
        A[U11, S21] = w * 2.0 * z

        A[U22, T22] = w * 1.0
        A[U22, L11] = w * zz
        A[U22, L33] = w * xx
        A[U22, L13] = w * -2.0 * xz
        A[U22, S12] = w * -2.0 * z
        A[U22, S32] = w * 2.0 * x

        A[U33, T33] = w * 1.0
        A[U33, L11] = w * yy
        A[U33, L22] = w * xx
        A[U33, L12] = w * -2.0 * xy
        A[U33, S23] = w * -2.0 * x
        A[U33, S13] = w * 2.0 * y

        A[U12, T12] = w * 1.0
        A[U12, L33] = w * -xy
        A[U12, L23] = w * xz
        A[U12, L13] = w * yz
        A[U12, L12] = w * -zz
        A[U12, S11] = w * -z
        A[U12, S22] = w * z
        A[U12, S31] = w * x
        A[U12, S32] = w * -y

        A[U13, T13] = w * 1.0
        A[U13, L22] = w * -xz
        A[U13, L23] = w * xy
        A[U13, L13] = w * -yy
        A[U13, L12] = w * yz
        A[U13, S11] = w * y
        A[U13, S33] = w * -y
        A[U13, S23] = w * z
        A[U13, S21] = w * -x

        A[U23, T23] = w * 1.0
        A[U23, L11] = w * -yz
        A[U23, L23] = w * -xx
        A[U23, L13] = w * xy
        A[U23, L12] = w * xz
        A[U23, S22] = w * -x
        A[U23, S33] = w * x
        A[U23, S12] = w * y
        A[U23, S13] = w * -z

        if args.biso:
            A[U11, i+21] = w
            A[U22, i+21] = w
            A[U33, i+21] = w

    # build B
    B = np.zeros((natoms * 6))
    for i in xrange(natoms):
        B[i * 6] = uu[i * 3, i * 3]
        B[i * 6 + 1] = uu[i * 3 + 1, i * 3 + 1]
        B[i * 6 + 2] = uu[i * 3 + 2, i * 3 + 2]
        B[i * 6 + 3] = uu[i * 3, i * 3 + 1]
        B[i * 6 + 4] = uu[i * 3, i * 3 + 2]
        B[i * 6 + 5] = uu[i * 3 + 1, i * 3 + 2]

    # back calculated U_tls and compare to values in B.
    UTLS = np.dot(A, xv)
    D = UTLS - B
    SSE = sum((UTLS - B)**2)
    MSE = SSE / len(B)
    # pearson correlation coefficinet
    r2 = (np.dot(B,UTLS))**2/(np.sum(B**2) * np.sum(UTLS**2))

    for i in range(len(B)):
        print " %12.5f  %12.5f" % (B[i], UTLS[i])

    if args.biso:
        # follow the procedure in Section 3.1.2 of Painter & Merritt,
        # 2006, to make all the Biso values positive:
        #
        #  Note: this is just a local calculation here, for printing to
        #  the ouptut file: it doesn't affect the TLS matrices or ADP's
        #  that get written to the output pdb file.
        #
        #  first, need to get the "pure" TLS motion by setting all
        #   columns of a beyond 21 to zero:
        for row in range(len(B)):
            for column in range(natoms):
                A[row,column+21] = 0
        UTLSp = np.dot(A,xv)

        print "#   i     Uobs      Utls        Uiso"

        # use label indexing to avoid confusion!
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
            S11, S22, S33, S12, S13, S23, S21, S31, S32 = \
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, \
            15, 16, 17, 18, 19, 20)

        Uisomin = 9999.
        for i in range(natoms):
            Utls = (UTLSp[6*i] + UTLSp[6*i+1] + UTLSp[6*i+2])/3.
            Uobs = (B[6*i] + B[6*i+1] + B[6*i+2])/3.
            Uiso = Uobs - Utls
            print "%5d  %8.3f %8.3f %8.3f" % ( i, Uobs, Utls, Uiso )
            if Uiso < Uisomin:
                Uisomin = Uiso
        print "Uisomin is %12.5f" % Uisomin

    return UTLS, SSE, MSE, r2


def center_of_reaction_shift(T0, L0, S0, origin):
    """Calculate new tensors based on the center for reaction.
    This method creates a dictionary of the calculations:

    T^: T tensor in the coordinate system of L
    L^: L tensor in the coordinate system of L
    S^: S tensor in the coordinate system of L

    COR: Center of Reaction

    T',S',L': T,L,S tensors in original coordinate system
              with the origin shifted to the center of reaction.

    It then returns T', S', L' and the COR for printing out to file.
    """

    # LSMALL is the smallest magnitude of L before it is considered 0.0
    LSMALL = 0.5 * (np.pi / 180.0)**2
    rdict = {}

    rdict["T'"] = T0.copy()
    rdict["L'"] = L0.copy()
    rdict["S'"] = S0.copy()

    rdict["rT'"] = T0.copy()

    rdict["L1_eigen_val"] = 0.0
    rdict["L2_eigen_val"] = 0.0
    rdict["L3_eigen_val"] = 0.0

    rdict["L1_rmsd"] = 0.0
    rdict["L2_rmsd"] = 0.0
    rdict["L3_rmsd"] = 0.0

    rdict["L1_eigen_vec"] = np.zeros(3, float)
    rdict["L2_eigen_vec"] = np.zeros(3, float)
    rdict["L3_eigen_vec"] = np.zeros(3, float)

    rdict["RHO"] = np.zeros(3, float)
    rdict["COR"] = origin

    rdict["L1_rho"] = np.zeros(3, float)
    rdict["L2_rho"] = np.zeros(3, float)
    rdict["L3_rho"] = np.zeros(3, float)

    rdict["L1_pitch"] = 0.0
    rdict["L2_pitch"] = 0.0
    rdict["L3_pitch"] = 0.0

    rdict["Tr1_eigen_val"] = 0.0
    rdict["Tr2_eigen_val"] = 0.0
    rdict["Tr3_eigen_val"] = 0.0

    rdict["Tr1_rmsd"] = 0.0
    rdict["Tr2_rmsd"] = 0.0
    rdict["Tr3_rmsd"] = 0.0

    # set the L tensor eigenvalues and eigenvectors
    (L_evals, RL) = linalg.eig(L0)
    RL = RL.T
    L1, L2, L3 = L_evals

    good_L_eigens = []

    if np.allclose(L1, 0.0) or isinstance(L1, complex):
        L1 = 0.0
    else:
        good_L_eigens.append(0)

    if np.allclose(L2, 0.0) or isinstance(L2, complex):
        L2 = 0.0
    else:
        good_L_eigens.append(1)

    if np.allclose(L3, 0.0) or isinstance(L3, complex):
        L3 = 0.0
    else:
        good_L_eigens.append(2)

    # no good L eigenvalues
    if len(good_L_eigens) == 0:
        return T0, L0, S0, origin

    # one good eigenvalue -- reconstruct RL about it
    elif len(good_L_eigens) == 1:
        i = good_L_eigens[0]
        evec = RL[i]

        RZt = np.transpose(AtomMath.rmatrixz(evec))
        xevec = np.dot(RZt, np.array([1.0, 0.0, 0.0], float))
        yevec = np.dot(RZt, np.array([0.0, 1.0, 0.0], float))

        if i == 0:
            RL[1] = xevec
            RL[2] = yevec
        elif i == 1:
            RL[0] = xevec
            RL[2] = yevec
        elif i == 2:
            RL[0] = xevec
            RL[1] = yevec

    # two good eigenvalues -- reconstruct RL about them
    elif len(good_L_eigens) == 2:
        i = good_L_eigens[0]
        j = good_L_eigens[1]

        xevec = AtomMath.normalize(np.cross(RL[i], RL[j]))
        for k in range(3):
            if k == i:
                continue
            if k == j:
                continue
            RL[k] = xevec
            break

    rdict["L1_eigen_val"] = L1
    rdict["L2_eigen_val"] = L2
    rdict["L3_eigen_val"] = L3

    if L1 > 0:
        rdict["L1_rmsd"] = np.sqrt(L1)
    else:
        rdict["L1_rmsd"] = 0
    if L2 > 0:
        rdict["L2_rmsd"] = np.sqrt(L2)
    else:
        rdict["L2_rmsd"] = 0
    if L3 > 0:
        rdict["L3_rmsd"] = np.sqrt(L3)
    else:
        rdict["L3_rmsd"] = 0

    rdict["L1_eigen_vec"] = RL[0].copy()
    rdict["L2_eigen_vec"] = RL[1].copy()
    rdict["L3_eigen_vec"] = RL[2].copy()

    # begin tensor transformations which depend upon
    # the eigenvectors of L0 being well-determined
    # make sure RLt is right-handed
    if np.allclose(linalg.det(RL), -1.0):
        I = np.identity(3, float)
        I[0, 0] = -1.0
        RL = np.dot(I, RL)

    if not np.allclose(linalg.det(RL), 1.0):
        return T0, L0, S0, origin

    RLt = np.transpose(RL)

    # carrot-L tensor (tensor WRT principal axes of L)
    cL = np.dot(np.dot(RL, L0), RLt)
    rdict["L^"] = cL.copy()

    # carrot-T tensor (T tensor WRT principal axes of L)
    cT = np.dot(np.dot(RL, T0), RLt)
    rdict["T^"] = cT.copy()

    # carrot-S tensor (S tensor WRT principal axes of L)
    cS = np.dot(np.dot(RL, S0), RLt)
    rdict["S^"] = cS.copy()

    # ^rho: the origin-shift vector in the coordinate system of L
    L23 = L2 + L3
    L13 = L1 + L3
    L12 = L1 + L2

    # shift for L1
    if not np.allclose(L1, 0.0) and abs(L23) > LSMALL:
        crho1 = (cS[1, 2] - cS[2, 1]) / L23
    else:
        crho1 = 0.0

    if not np.allclose(L2, 0.0) and abs(L13) > LSMALL:
        crho2 = (cS[2, 0] - cS[0, 2]) / L13
    else:
        crho2 = 0.0

    if not np.allclose(L3, 0.0) and abs(L12) > LSMALL:
        crho3 = (cS[0, 1] - cS[1, 0]) / L12
    else:
        crho3 = 0.0

    crho = np.array([crho1, crho2, crho3], float)
    rdict["RHO^"] = crho.copy()

    # rho: the origin-shift vector in orthogonal coordinates
    rho = np.dot(RLt, crho)
    rdict["RHO"] = rho
    rdict["COR"] = origin + rho

    # set up the origin shift matrix PRHO WRT orthogonal axes
    PRHO = np.array([[0.0, rho[2], -rho[1]],
                     [-rho[2], 0.0, rho[0]],
                     [rho[1], -rho[0], 0.0]], float)

    # set up the origin shift matrix cPRHO WRT libration axes
    cPRHO = np.array([[0.0, crho[2], -crho[1]],
                      [-crho[2], 0.0, crho[0]],
                      [crho[1], -crho[0], 0.0]], float)

    # calculate transpose of cPRHO, ans cS
    cSt = np.transpose(cS)
    cPRHOt = np.transpose(cPRHO)

    rdict["L'^"] = cL.copy()

    # calculate S'^ = S^ + L^*pRHOt
    cSp = cS + np.dot(cL, cPRHOt)
    rdict["S'^"] = cSp.copy()

    # calculate T'^ = cT + cPRHO*S^ + cSt*cPRHOt + cPRHO*cL*cPRHOt *
    cTp = cT + np.dot(cPRHO, cS) + np.dot(cSt, cPRHOt) +\
        np.dot(np.dot(cPRHO, cL), cPRHOt)
    rdict["T'^"] = cTp.copy()

    # transpose of PRHO and S
    PRHOt = np.transpose(PRHO)
    St = np.transpose(S0)

    # calculate S' = S + L*PRHOt
    Sp = S0 + np.dot(L0, PRHOt)
    rdict["S'"] = Sp

    # calculate T' = T + PRHO*S + St*PRHOT + PRHO*L*PRHOt
    Tp = T0 + np.dot(PRHO, S0) + np.dot(St, PRHOt) +\
        np.dot(np.dot(PRHO, L0), PRHOt)
    rdict["T'"] = Tp

    # now calculate the TLS motion description using 3 non
    # intersecting screw axes, with one

    # libration axis 1 shift in the L coordinate system
    if abs(L1) > LSMALL:
        cL1rho = np.array([0.0, -cSp[0, 2] / L1, cSp[0, 1] / L1], float)
    else:
        cL1rho = np.zeros(3, float)

    # libration axis 2 shift in the L coordinate system
    if abs(L2) > LSMALL:
        cL2rho = np.array([cSp[1, 2] / L2, 0.0, -cSp[1, 0] / L2], float)
    else:
        cL2rho = np.zeros(3, float)

    # libration axis 2 shift in the L coordinate system
    if abs(L3) > LSMALL:
        cL3rho = np.array([-cSp[2, 1] / L3, cSp[2, 0] / L3, 0.0], float)
    else:
        cL3rho = np.zeros(3, float)

    # libration axes shifts in the original orthogonal coordinate system
    rdict["L1_rho"] = np.dot(RLt, cL1rho)
    rdict["L2_rho"] = np.dot(RLt, cL2rho)
    rdict["L3_rho"] = np.dot(RLt, cL3rho)

    # calculate screw pitches (A*R / R*R) = (A/R)
    if abs(L1) > LSMALL:
        rdict["L1_pitch"] = cS[0, 0] / L1
    else:
        rdict["L1_pitch"] = 0.0

    if L2 > LSMALL:
        rdict["L2_pitch"] = cS[1, 1] / L2
    else:
        rdict["L2_pitch"] = 0.0

    if L3 > LSMALL:
        rdict["L3_pitch"] = cS[2, 2] / L3
    else:
        rdict["L3_pitch"] = 0.0

    # now calculate the reduction in T for the screw rotation axes
    cTred = cT.copy()

    for i in (0, 1, 2):
        for k in (0, 1, 2):
            if i == k:
                continue
            if abs(cL[k, k]) > LSMALL:
                cTred[i, i] -= (cS[k, i]**2) / cL[k, k]

    for i in (0, 1, 2):
        for j in (0, 1, 2):
            for k in (0, 1, 2):
                if j == i:
                    continue
                if abs(cL[k, k]) > LSMALL:
                    cTred[i, j] -= (cS[k, i] * cS[k, j]) / cL[k, k]

    # rotate the newly calculated reduced-T tensor from the carrot
    # coordinate system (coordinate system of L) back to the structure
    # coordinate system
    Tr = np.dot(np.dot(RLt, cTred), RL)
    rdict["rT'"] = Tr

    Tr1, Tr2, Tr3 = linalg.eigvals(Tr)

    if np.allclose(Tr1, 0.0) or isinstance(Tr1, complex):
        Tr1 = 0.0
    if np.allclose(Tr2, 0.0) or isinstance(Tr2, complex):
        Tr2 = 0.0
    if np.allclose(Tr3, 0.0) or isinstance(Tr3, complex):
        Tr3 = 0.0

    rdict["Tr1_eigen_val"] = Tr1
    rdict["Tr2_eigen_val"] = Tr2
    rdict["Tr3_eigen_val"] = Tr3

    if Tr1 > 0:
        rdict["Tr1_rmsd"] = np.sqrt(Tr1)
    else:
        rdict["Tr1_rmsd"] = 0
    if Tr2 > 0:
        rdict["Tr2_rmsd"] = np.sqrt(Tr2)
    else:
        rdict["Tr2_rmsd"] = 0
    if Tr3 > 0:
        rdict["Tr3_rmsd"] = np.sqrt(Tr3)
    else:
        rdict["Tr3_rmsd"] = 0

    T = rdict["T'"].copy()
    L = rdict["L'"].copy()
    S = rdict["S'"].copy()
    com = rdict["COR"].copy()

    return T, L, S, com


def run(infile, prefix, covariance_file):
    # read in pdb file (parmed)
    pdb = read_PDB(infile)
    natoms = len(pdb.atoms)
    assert pdb.box is not None, "ERROR: CRYST1 record required in PDB file!"
    # get centroid
    com = center_of_mass(pdb)

    if not covariance_file:
        # set A
        A = set_A(pdb, com)
        # set B
        B = set_B(pdb)
        # solve Ax=B
        x = solve_SVD(A, B)
        # calc new ADP and residual and print residual
        print "# b, A.x for individual atom covariance fit:"
        UTLS, SSE, MSE, r2 = calc_ADP_from_TLS(x, A, B, natoms)

    else:
        # expanded A with 9*n(n+1)/2-3n equations
        A_cov = set_A_cov(pdb, com)
        # 9*n(n+1)/2-3n covariances
        B_cov, uu = set_B_cov(pdb, covariance_file)
        # 21 TLS parameters
        x = solve_SVD(A_cov, B_cov)

        # get the statistics of the fit to the full covariance matrix:
        print "# b, A.x for full covariance matrix fit:"
        calc_ADP_from_TLS(x, A_cov, B_cov, natoms)

        # calc new ADP and residual and print residual
        print "# b, A.x for individual atoms, extracted from full covariance matrix fit:"
        UTLS, SSE, MSE, r2 = calc_ADP_from_TLS_cov(pdb, com, uu, x)

    # print statistics
    print "#      *** ADP FIT STATISTICS FOR INDIVIDUAL ATOMS ***"
    print "# The ADP Sum of Squared Errors is %6.4f." % SSE
    print "# The ADP Mean Squared Error is %6.4f." % MSE
    print "# The ADP Pearson R-square is %6.4f." % r2

    # populate TLS matrices
    T, L, S = set_TLS(x)
    cor = com

    # optionally, shift TLS to center of reaction
    if not args.nocor:
        print "# Shifting TLS to the center of reaction"
        T, L, S, cor = center_of_reaction_shift(T, L, S, com)

    # write TLS to file
    write_TLS(T, L, S, cor, pdb, prefix)
    write_TLS_remark3(T, L, S, cor, pdb, prefix)

    # write out pdb file with new ADP/Bfactors
    for i, atom in enumerate(pdb.atoms):
        atom.anisou = UTLS[i * 6:i * 6 + 6]
        Bfac, Ueq = Adp2B(np.array(pdb.box), atom.anisou, method=2)
        atom.bfactor = Bfac
    write_PDB(pdb, "%s.pdb1" % prefix, write_anisou=True)

    # put output pdb file and REMARK 3 cards into the same file
    os.system( 'cat %s.remark3.dat %s.pdb1 > %s.pdb' % ( prefix, prefix,
           prefix ) )
    os.system( '/bin/rm -f %s.remark3.dat %s.pdb1' % (prefix, prefix))

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="input PDB file (min 4 atoms).")
    parser.add_argument(
        "-o",
        "--prefix",
        help="prefix for output .pdb and .dat files",
        default="out")
    parser.add_argument(
        "-cov",
        "--covariance",
        help="covariance matrix from cpptraj (optional input)",
        default=None)
    parser.add_argument(
        "-nocor",
        "--no-center-of-reaction-shift", dest='nocor', 
        default=False, action='store_true',
        help="skip moving to the center of reaction")
    parser.add_argument(
        "-biso",
        "--biso", dest='biso', 
        default=False, action='store_true',
        help="fit TLS + isoB model")
    args = parser.parse_args()
    if not args.infile:
        parser.print_help()
        sys.exit(1)
    assert os.path.isfile(args.infile), "ERROR: input PDB file not found!"
    if args.covariance:
        assert os.path.isfile(
            args.covariance), "ERROR: covariance matrix file file not found!"
    print "#  Input pdb file:  %s" % args.infile
    print "#  Output prefix :  %s" % args.prefix
    if args.covariance:
        print "#  Input cov file:  %s" % args.covariance
   
    run(args.infile, args.prefix, args.covariance)

