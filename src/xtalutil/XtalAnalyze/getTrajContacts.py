
########################################################################
# SCRIPT TO GET CRYSTAL CONTACTS FOR A SUPERCELL                       #
# FOR USAGE TYPE -H ON THE COMMANDLINE:                                #
#      phenix.python getTrajContacts.py -h                             #
########################################################################

from libtbx import group_args
import sys
import os
import numpy as np
from iotbx import file_reader


########################################################################
# This is the basic phenix crystal contact finder that I got from James#
# Fraser. Returns a list of contacts of each atom with all other atoms #
# within a certain cut-off.                                            #
########################################################################
def find_crystal_contacts(xray_structure,
                          pdb_atoms,  # atom_with_labels, not atom!
                          selected_atoms=None,
                          distance_cutoff=3.5,
                          ignore_same_asu=True,
                          ignore_waters=True):
    from scitbx.array_family import flex
    sites_frac = xray_structure.sites_frac()
    unit_cell = xray_structure.unit_cell()
    pair_asu_table = xray_structure.pair_asu_table(
        distance_cutoff=distance_cutoff)
    pair_sym_table = pair_asu_table.extract_pair_sym_table()
    contacts = []
    if (selected_atoms is None):
        selected_atoms = flex.bool(len(pdb_atoms), True)
    for i_seq, pair_sym_dict in enumerate(pair_sym_table):
        if (not selected_atoms[i_seq]):
            continue
        site_i = sites_frac[i_seq]
        atom_i = pdb_atoms[i_seq]
        resname_i = atom_i.resname
        atmname_i = atom_i.name
        chainid_i = atom_i.chain_id
        for j_seq, sym_ops in list(pair_sym_dict.items()):
            site_j = sites_frac[j_seq]
            atom_j = pdb_atoms[j_seq]
            resname_j = atom_j.resname
            atmname_j = atom_j.name
            chainid_j = atom_j.chain_id
            for sym_op in sym_ops:
                if sym_op.is_unit_mx():
                    if ignore_same_asu:
                        continue
                    #~ elif (chainid_i == chainid_j) :
                        #~ continue
                if (resname_j in ["HOH", "WAT"] and ignore_waters):
                    continue
                site_ji = sym_op * site_j
                distance = unit_cell.distance(site_i, site_ji)
                contacts.append((i_seq, j_seq, sym_op, distance))
                # print resname_i, atmname_i, resname_j, atmname_j, str(sym_op), distance
    return contacts

########################################################################
#                                                                      #
# This takes the contacts list from find_crystal_contacts() which has  #
# contact entries for each atom and returns all_residues which is a    #
# nested list that has entries for each residue.                       #
# Returns:                                                             #
#  all_residues[i] - each entry corresponds to one residue             #
#  all_residues[i][0] - list of the residue name, number...            #
#  all_residues[i][1] - list of all atoms interacting with the residue #
#  all_residues[i][1][i] - list of atom no., symmetry op., ...         #
#                                                                      #
########################################################################


def find_crystal_contacts_by_residue(xray_structure,
                                     pdb_hierarchy,
                                     **kwds):
    contacts_by_residue = {}
    atoms = list(pdb_hierarchy.atoms_with_labels())
    contacts = find_crystal_contacts(xray_structure, atoms, **kwds)
    for (i_seq, j_seq, sym_op, distance) in contacts:
        atom_rec = atoms[i_seq].fetch_labels()
        residue_key = (atom_rec.chain_id, atom_rec.resname, atom_rec.resid(),
                       atom_rec.altloc)
        if (not residue_key in contacts_by_residue):
            contacts_by_residue[residue_key] = []
        contacts_by_residue[residue_key].append((j_seq, sym_op, distance))
    all_residues = []
    for chain in pdb_hierarchy.models()[0].chains():
        chain_id = chain.id
        for residue_group in chain.residue_groups():
            resid = residue_group.resid()
            for atom_group in residue_group.atom_groups():
                resname = atom_group.resname
                altloc = atom_group.altloc
                residue_key = (chain_id, resname, resid, altloc)
                residue_contacts = contacts_by_residue.get(residue_key, [])
                all_residues.append((residue_key, residue_contacts))

    return all_residues

########################################################################
# For residue id in a super-cell, get residue id in original asymmetric#
# unit, symmetry operation, and unit cell translation.                 #
# Arguments:                                                           #
#     prop_vec: 3x1 array of the x,y,z propagations used to create the #
#               the super-cell                                         #
#     nres: int number of residues in asymmetric unit                  #
#     nsymop: int number of symmetry operations in the space group     #
#     resid: int query residue number in the super-cell                #
# Returns:                                                             #
#     sc_trans: 3x1 tuple of x,y,z propagations to arrive at the query #
#               residue's unit cell within the super-cell              #
#     symop: int number of the symop used to get to the query residue  #
#            (starts with 0)                                           #
#     asymresid: int corresponding residue id in the asymmetric unit   #
########################################################################


def SC_mapping(prop_vec, nres, nsymop, resid):
    from cctbx_sgtbx_ext import rt_mx
    xnres = prop_vec[2] * prop_vec[1] * nres * nsymop
    ynres = prop_vec[2] * nres * nsymop
    znres = nres * nsymop
    (xmoves, resid) = divmod(resid - 1, xnres)
    (ymoves, resid) = divmod(resid, ynres)
    (zmoves, resid) = divmod(resid, znres)
    (sym_op, resid) = divmod(resid, nres)
    tr_op = rt_mx('x+%d,y+%d,z+%d' % (xmoves, ymoves, zmoves))
    return tr_op, sym_op, resid + 1


def is_same_asu(nres, resid1, resid2):
    if divmod(resid1 - 1, nres)[0] == divmod(resid2 - 1, nres)[0]:
        return True
    else:
        return False

########################################################################
#                                                                      #
# Takes the all_residues list, converts the symmetry operations        #
# relative to the supercell as P1 unit cell to symmetry operations     #
# relative to the ASU and the true space group symmetry. Returns       #
# supercell_contacts                                                   #
########################################################################


def find_supercell_contacts_by_residue(residue_contacts, atoms, prop_vec,
                                       nres, nsymop, rt_mx_matrices, ignore_waters=True):
    supercell_contacts = {}
    # skip all solvent ions added that does not belong to any asymmetric unit
    # by calculating maxresidue and skipping residues above that
    maxres = nres * prop_vec[0] * prop_vec[1] * prop_vec[2] * nsymop
    # loop over each residue in residue_contacts
    for (residue, contacts) in residue_contacts:
        i_resid = int(residue[2])
        i_resname = residue[1]
        if (ignore_waters and i_resname in ["HOH", "WAT"]):
            continue
        if i_resid >= maxres:
            continue

        # residue_contacts may have multiple entries between two residues due to
        # different atoms of the same residues. Duplicate_test{} is used
        # to filter this so that only one contact between each pair of residues
        # is returned with the shortest atom to atom distance reported.
        # sc_sym_op is the translation that relates the original supercell to
        # the supercell where residue_j resides when the contact occurs (remember
        # that the supercell is being treated as a P1 'unit cell'
        duplicate_test = {}
        for contact in contacts:
            sc_sym_op = contact[1]
            j_seq = contact[0]
            j_resid = int(atoms[j_seq].fetch_labels().resid())
            j_resname = atoms[j_seq].fetch_labels().resname
            dist = contact[2]
            if j_resid > maxres:
                continue

            # unit_mx sc_sym_op means the contact is in the original supercell
            # if same asu in same supercell, it's an "intra_asu" contact so
            # ignore
            if is_same_asu(nres, i_resid, j_resid):
                if sc_sym_op.is_unit_mx():
                    continue
            if (ignore_waters and j_resname in ["HOH", "WAT"]):
                continue

            # populate duplicate_test{} with only one instance of a given
            # residue pair, keeping the shortest distance value
            if (j_resid, j_resname, sc_sym_op) in list(duplicate_test.keys()):
                if dist < duplicate_test[(j_resid, j_resname, sc_sym_op)]:
                    duplicate_test[(j_resid, j_resname, sc_sym_op)] = dist
            else:
                duplicate_test[(j_resid, j_resname, sc_sym_op)] = dist

        # for contact residue pair, get residue's number in original asu,
        # symop number to arrive at the residue's asu (this is the ordinal number
        # of the symop from the rt_mx_matrices object (SYMM cards in PDB file)
        # applied to the resiude in the original ASU to bring it to it's current
        # ASU), and the translation operation to bring it to its current unit cell
        for (j_resid, j_resname, sc_sym_op) in list(duplicate_test.keys()):
            dist = duplicate_test[(j_resid, j_resname, sc_sym_op)]
            i_transop, symop1, i_asymresid = SC_mapping(prop_vec, nres, nsymop, i_resid)
            j_transop, symop2, j_asymresid = SC_mapping(prop_vec, nres, nsymop, j_resid)

            # both residues within the original supercell, so get the translation
            # that relates j's unit_cell to i's unit_cell
            if sc_sym_op.is_unit_mx():
                ij_transop = j_transop.multiply(i_transop.inverse())
            # res_j outside the original supercell
            else:
                # get the translation operation of the supercell where the contact
                # residue would reside relative to the original supercell
                t = sc_sym_op.t().as_double()
                # augment by the translations inherent in the supercell itself
                t = [t[0] * prop_vec[0], t[1] * prop_vec[1], t[2] * prop_vec[2]]
                # convert back to rt_mx
                t = [int(i) for i in t]
                from cctbx_sgtbx_ext import rt_mx, tr_vec

                # t is the translation vector to go from an asu in the original
                #           supercell to an asu in the supercell where the
                #           contact would occur
                # j_transop is the translation from the original unit cell to the
                #           unit cell in the supercell where resid_j is located
                # to get the translation necessary to go from resid_i unit cell to
                #           resid_j supercell we do:
                #
                #                 t+j_transop - i_transop
                #
                t = rt_mx(tr_vec(t, tr_den=1))
                j_transop = j_transop.multiply(t)
                ij_transop = j_transop.multiply(i_transop.inverse())


# THIS SHOULD WORK BUT DOESNT
            # get the symmetry operation that relates the asu of resid_j to
            # the asu of resid_i:
            ##
            # Rij*Ri=Rj
            # Rij=Rj*Ri^-1
            ##
            # i_r=rt_mx_matrices[symop1]
            # j_r=rt_mx_matrices[symop2]
            # ij_r=j_r.multiply(i_r.inverse())

            # now add the translation that relates unit cell of resiude_j to
            # unit cell of residue_i:
            ##
            # sym_op=Rij+transop
            ##
            # sym_op=ij_transop.multiply(ij_r)

            # this symop is relative to the asu of residue_i. Apply the inverse
            # of the symmetry operation that creates residue_i's asu to obtain
            # the symop relative to the original asu.
            # transop=i_r.inverse().multiply(sym_op)
# THIS SHOULD WORK BUT DOESNT

# THIS WORKS ON MANY SPACEGROUPS BUT NOT P212121
            # The transop translation is relative to the asu of resid_i. Rotate
            # it by the inverse of the rotation part of the rot/trans matrix
            # used to create the asu of resid_i to obtain the translation relative
            # to the original asu.
            rot = rt_mx_matrices[symop1]
            transop = rot.inverse().r().multiply(ij_transop.t())
            transop = rt_mx(transop)

            # Obtain the rot/trans matrix that relates the the asu of resid_j to
            # the asu of resid_i.
            i_r = rt_mx_matrices[symop1]
            j_r = rt_mx_matrices[symop2]
            ij_r = j_r.multiply(i_r.inverse())

            # Theoretically this rotation should now be rotated by the inverse of
            # the rotation used to get to asu of resid_i to get the rotation
            # relative to the original ASU. But this does not work; works without
            # it. Don't know why.
            # TEST
            #~ ij_r=i_r.inverse().multiply(ij_r)

            # Add the translation relative to the original ASU to the symmetry op
            # relative to the original ASU
            transop = transop.multiply(ij_r)
# THIS WORKS ON MANY SPACEGROUPS BUT NOT P212121

            # populate supercell_contacts{}. To avoid duplicates chose by highest
            # x trans op, then y then z
            if transop.as_xyz() == transop.inverse().as_xyz():
                if i_asymresid <= j_asymresid:
                    i_key = (
                        i_asymresid,
                        i_resname,
                        j_asymresid,
                        j_resname,
                        transop.as_xyz())
                else:
                    i_key = (
                        j_asymresid,
                        j_resname,
                        i_asymresid,
                        i_resname,
                        transop.as_xyz())
            elif transop.t().as_double()[0] > transop.inverse().t().as_double()[0]:
                i_key = (i_asymresid, i_resname, j_asymresid, j_resname, transop.as_xyz())
            elif transop.t().as_double()[0] < transop.inverse().t().as_double()[0]:
                i_key = (
                    j_asymresid,
                    j_resname,
                    i_asymresid,
                    i_resname,
                    transop.inverse().as_xyz())
            else:
                if transop.t().as_double()[1] > transop.inverse().t().as_double()[1]:
                    i_key = (
                        i_asymresid,
                        i_resname,
                        j_asymresid,
                        j_resname,
                        transop.as_xyz())
                elif transop.t().as_double()[1] < transop.inverse().t().as_double()[1]:
                    i_key = (
                        j_asymresid,
                        j_resname,
                        i_asymresid,
                        i_resname,
                        transop.inverse().as_xyz())
                else:
                    if transop.t().as_double()[2] > transop.inverse().t().as_double()[2]:
                        i_key = (
                            i_asymresid,
                            i_resname,
                            j_asymresid,
                            j_resname,
                            transop.as_xyz())
                    elif transop.t().as_double()[2] < transop.inverse().t().as_double()[2]:
                        i_key = (
                            j_asymresid,
                            j_resname,
                            i_asymresid,
                            i_resname,
                            transop.inverse().as_xyz())
                    else:
                        print("Warning: maybe something is wrong. Transop is: ")
                        print(transop.as_xyz())
                        import code
                        code.interact(local=locals())
                        sys.exit()

            #~ print i_resid, j_resid, i_asymresid, j_asymresid, transop

            # The residue numbers, names and symmetry operation uniquely identify
            # each contact. Populate supercell_contacts dictionary object. Keys
            # are the cotnact identified. Values are list of distances.
            if i_key in list(supercell_contacts.keys()):
                supercell_contacts[i_key].append(dist)
            else:
                supercell_contacts[i_key] = [dist]
    return supercell_contacts

########################################################################
# Report                                                               #
########################################################################


def report_supercell_contacts(supercell_contacts):
    transops = set([i[4] for i in list(supercell_contacts.keys())])
    for transop in transops:
        tmp_l = [(k, v) for k, v in supercell_contacts.items() if k[4] == transop]
        print("S%12s %3d|\n" % (transop, len(tmp_l)), end=' ')
        for h, i in enumerate(tmp_l):
            print("                  %4d %3s %4d %3s %3d" \
                % (i[0][0], i[0][1], i[0][2], i[0][3], len(i[1])))


########################################################################
# This returns a list of rt_mx objects containing the rot/tr operations#
# to create the unit cell from the asu for a given spacegroup          #
# Input:                                                               #
#    sg - string containing the name of the space group (eg. 'P1211')  #
# Returns:                                                             #
#    rt_mx_matrices - list of rt_mx objects                            #
########################################################################
#~ def get_symm(sg):
    #~ from cctbx.sgtbx import space_group_info
    #~ i=space_group_info(sg)
    #~ rt_mx_matrices=[rt_mx for rt_mx in i.group().smx()]
    #~ return rt_mx_matrices

def get_symm(sg):
    from cctbx_sgtbx_ext import rt_mx
    if sg == 'P212121':
        rot0 = rt_mx("x,y,z")
        rot1 = rt_mx("-x+1/2,-y,z+1/2")
        rot2 = rt_mx("-x,y+1/2,-z+1/2")
        rot3 = rt_mx("x+1/2,-y+1/2,-z")
        rt_mx_matrices = (rot0, rot1, rot2, rot3)

    elif sg == 'P1211':
        rot0 = rt_mx("x,y,z")
        rot1 = rt_mx("-x,y+1/2,-z")
        rt_mx_matrices = (rot0, rot1)

    elif sg == 'P121':
        rot0 = rt_mx("x,y,z")
        rot1 = rt_mx("-x,y,-z")
        rt_mx_matrices = (rot0, rot1)

    elif sg == 'P222':
        rot0 = rt_mx("x,y,z")
        rot1 = rt_mx("-x,-y,z")
        rot2 = rt_mx("-x,y,-z")
        rot3 = rt_mx("x,-y,-z")
        rt_mx_matrices = (rot0, rot1, rot2, rot3)

    else:
        print("%s not found\n" % sg)
        sys.exit()

    return rt_mx_matrices


if __name__ == "__main__":
    import argparse
    msg = """Eg. usage:
             phenix.python getTrajContacts.py 1d23X.pdb -n 20 -sg P212121 -c 2.5 -px 1 -py 1 -pz 1
             """
    #~ def msg(name=None):
    #~ return '''program.py
    #~ [-a, Pass argument a]
    #~ [-b, Pass argument b]
    #~ [-c, Pass argument c]
    #~ [-d, Pass argument d]
    #~ comment
    #~ more comment
    #~ '''

    parser = argparse.ArgumentParser(
        description='''This script takes a pdb supercell structure that consists
    of a unit cell propagated in x,y,z directions and returns for each
    crystal interface of the asymmetric unit, the number of unique contacts
    in that interface. Installation of cctbx or use of phenix.python is
    required. Note that the supercell PDB file must be provided with a
    CRYST1 record describing the supercell as a P1 unit cell.''',
        epilog=msg)

    parser.add_argument("file", help="PDB file. Must include CRYST1 line"
                        "with the provided supercell treated as a P1 unit"
                        "cell.")
    parser.add_argument("-n", help="No. of residues in asymme"
                        "tric unit.", type=int, required=True)
    parser.add_argument("-sg", help="Space group name, no"
                        " spaces. Eg. \'P1211\'", required=True)
    parser.add_argument("-c", help="Angstrom cutoff for "
                        "contacts. Default=3.", type=float, default=3.)
    parser.add_argument("-px", help="unit cell propagations along a-axis "
                        "to form supercell", type=int, action='append', dest='prop', required=True)
    parser.add_argument("-py", help="unit cell propagations along a-axis "
                        "to form supercell", type=int, action='append', dest='prop', required=True)
    parser.add_argument("-pz", help="unit cell propagations along a-axis "
                        "to form supercell", type=int, action='append', dest='prop', required=True)
    args = parser.parse_args()

    pdb_file = args.file  # CRYST1 record w/supercell box required
    nres = args.n  # n residues in asymmetric unit
    sg = args.sg  # space group
    cutoff = args.c  # angstrom cutoff to find contacts
    prop_vec = args.prop  # propagation forming supercell

    pdb_in = file_reader.any_file(pdb_file).file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    xrs = pdb_in.xray_structure_simple()
    rt_mx_matrices = get_symm(sg)
    residue_contacts = find_crystal_contacts_by_residue(xrs, pdb_hierarchy,
                                                        distance_cutoff=cutoff, ignore_same_asu=False, ignore_waters=False)
    atoms = list(pdb_hierarchy.atoms_with_labels())
    supercell_contacts = find_supercell_contacts_by_residue(residue_contacts,
                                                            atoms, prop_vec, nres, len(rt_mx_matrices), rt_mx_matrices, ignore_waters=True)
    report_supercell_contacts(supercell_contacts)


# THE END

#======================================================================#
# 1rpg nres=124, nsymop2, P1211
# 1d23 nres=20, nsymop4, P212121
# 1jmt nres=124, P222
# 3pro nres= , P1211
# 4acu nres= , P212121

    #~ prop_vec=[2,2,2]
    #~ nres=20             #n residues in asymmetric unit
    #~ sg='P212121'            #space group
    #~ cutoff=3           #angstrom cutoff to find contacts
    #~ pdb_file = sys.argv[1]    #CRYST1 record w/supercell box required

    #~ import code; code.interact(local=locals())

    #~ l=sorted(supercell_contacts.iteritems(), key=lambda x: x[0][4] )
    #~ for i in l:
    #~ print i

    #~ if i_resname=="ASN" and i_resid==324:
    #~ print duplicate_test
    #~ for key in duplicate_test.keys():
    #~ print key[2]
    #~ debug=True

    #~ from cctbx_sgtbx_ext import rt_mx
    #~ rot=rot between symop1 and symop2

    #~ #P1211
    #~ rot0=rt_mx("x,y,z")
    #~ rot1=rt_mx("-x,y+1/2,-z")
    #~ rt_mx_matrices=(rot0,rot1)
    #~
    #~ for i in symm:
    #~ print i

    #~ #P121
    #~ rot0=rt_mx("x,y,z")
    #~ rot1=rt_mx("-x,y,-z")
    #~ rt_mx_matrices=(rot0,rot1)

    #~ #P212121
    #~ rot0=rt_mx("x,y,z")
    #~ rot1=rt_mx("-x+1/2,-y,z+1/2")
    #~ rot2=rt_mx("-x,y+1/2,-z+1/2")
    #~ rot3=rt_mx("x+1/2,-y+1/2,-z")
    #~ rt_mx_matrices=(rot0,rot1,rot2,rot3)

    #~ #P222
    #~ rot0=rt_mx("x,y,z")
    #~ rot1=rt_mx("-x,-y,z")
    #~ rot2=rt_mx("-x,y,-z")
    #~ rot3=rt_mx("x,-y,-z")
    #~ rt_mx_matrices=(rot0,rot1,rot2,rot3)

    #~ if i_resid==139:
    #~ import code; code.interact(local=dict(globals(), **locals()))

    # if symop2-symop1==1:
    # ~ transop=rot.inverse().multiply(transop.inverse())
    #transop = rot1.multiply(transop)
    # elif symop2-symop1==-1:
    # ~ rot=rot.inverse()
    # transop=rot1.inverse().multiply(transop)
    # ~ transop = transop.multiply(rot)
    # elif symop1-symop2==0:
    # pass

    #~ if i_asymresid in [124,32]:
    #~ print i_resid, j_resid, i_asymresid, j_asymresid, transop

    #~ if transop.as_xyz()=="-x-1,y+1/2,-z":
    #~ print i_resid, j_resid, i_asymresid, j_asymresid
    #~ if transop.inverse().as_xyz()=="-x-1,y+1/2,-z":
    #~ print i_resid, j_resid, i_asymresid, j_asymresid
    #~ if i_asymresid in [34, 66]:
    #~ print i_resid, j_resid, i_asymresid, j_asymresid, transop

    # This was the old way. It would avoid duplicates in the same frame but
    # not between frames (eg. would have (x+1,y,z) iface in one frame
    # and (x-1,y,z) contacts in another.
    #~ if transop.inverse().as_xyz() in [i[4] for i in supercell_contacts.keys()]:
    #~ i_key=(j_asymresid, j_resname, i_asymresid, i_resname,transop.inverse().as_xyz())
    #~ else:
    #~ i_key=(i_asymresid, i_resname, j_asymresid, j_resname,transop.as_xyz())

    #~ import code; code.interact(local=dict(globals(), **locals()))

    #~ for i in rt_mx_matrices:
    #~ print i
