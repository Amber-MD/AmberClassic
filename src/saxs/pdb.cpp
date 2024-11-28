#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include "saxsDS.h"
#include "pdb.h"
#include "const.h"

/////////////////
// Read pdb file
/////////////////
void read_pdb (const std::string &pdb_file,
               std::vector<coordinate> &pdb_coord) {

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        while (getline(PDBFILE, line))
            if ((line.find("ATOM") == 0) or (line.find("HETATM") == 0)) {
                coordinate coord;
                coord.x = atof (line.substr(30,8).c_str());
                coord.y = atof (line.substr(38,8).c_str());
                coord.z = atof (line.substr(46,8).c_str());
                coord.B_factor = atof (line.substr(60,6).c_str());

                std::string atm = line.substr(12,4);
                if (atm.substr(0,1) == " ")
                    atm = atm.substr(1,3);
                std::string type = atm.substr(0,1);
                if ((atm.substr(0,2) == "FE") or (atm.substr(0,2) == "Fe")) {
                    type = "Fe";
                    coord.r = 1.3;
                }

                if (type == "H")
                    coord.r = 1.2;
                else if (type == "C")
                    coord.r = 1.7;
                else if (type == "N")
                    coord.r = 1.55;
                else if (type == "O")
                    coord.r = 1.52;
                else if ((type == "S") or (type == "P"))
                    coord.r = 1.8;

                coord.type = type;
                coord.nHyd = 0;
                pdb_coord.push_back (coord);
            }
        PDBFILE.close();
    } else {
        std::cerr << "Unable to open file " << pdb_file << std::endl;
        exit (0);
    }
}

//////////////////////////////////////////////////
// Read trajectory file, each model to one vector
/////////////////////////////////////////////////
void read_pdb_weight (const std::string &pdb_file,
                      std::vector< std::vector <coordinate> > &snapshot,
                      std::vector<unsigned> &weight) {

    std::ifstream PDBFILE (pdb_file.c_str());
    if (PDBFILE.is_open()) {
        std::string line;
        std::vector<coordinate> model;
        while (getline(PDBFILE, line))
            if (line.find("MODEL") == 0)
                model.clear();          // Remove the previous model, ready to read the new model
            else if ((line.find("ATOM") == 0) or (line.find("HETATM") == 0)) {
                coordinate coord;
                coord.x = atof (line.substr(30,8).c_str());
                coord.y = atof (line.substr(38,8).c_str());
                coord.z = atof (line.substr(46,8).c_str());

                std::string atm = line.substr(12,4);
                if (atm.substr(0,1) == " ")
                    atm = atm.substr(1,3);

                std::string type;
                std::string ion = atm.substr(0,2);
                if ((ion == "Na") or (ion == "NA"))
                    type = "Na+";
                else if (ion == "K ")
                    type = "K+";
                else if ((ion == "Rb") or (ion == "RB"))
                    type = "Rb+";
                else if ((ion == "Cs") or (ion == "CS"))
                    type = "Cs+";
                else if ((ion == "Cl") or (ion == "CL"))
                    type = "Cl-";
                else if ((ion == "Br") or (ion == "BR"))
                    type = "Br-";
                else if ((ion == "Mg") or (ion == "MG"))
                    type = "Mg2+";
                else if (ion == "Ca")
                    type = "Ca2+";
                else if ((ion == "Sr") or (ion == "SR"))
                    type = "Sr2+";
                else if ((ion == "Ba") or (ion == "BA"))
                    type = "Ba2+";
                else {
                    type = atm.substr(0,1);
                    if ((line.find("WAT") != std::string::npos) or (line.find("SPC") != std::string::npos) or (line.find("T3P") != std::string::npos) or (line.find("T4E") != std::string::npos))
                        type += "w";        // Append "w" to differentiate water O and H atoms
                }
                coord.type = type;
                coord.nHyd = 0;
                model.push_back (coord);
            } else if ((line.find("ENDMDL") == 0) or (line.find("END") == 0))
                snapshot.push_back (model);
            else if (line.find("WEIGHT") == 0) {
                unsigned w = atoi (line.substr(9,6).c_str());
                // Assign weight to the current model
                weight.resize(snapshot.size()+1);
                weight.back() = w;
            }
        // In case there is no MODEL and ENDMDL keyword
        if (snapshot.size() == 0)
            snapshot.push_back (model);

        PDBFILE.close();
    } else {
        std::cerr << "Unable to open file " << pdb_file << std::endl;
        exit (0);
    }
    // Models must have the same size
    for (size_t i = 0; i < snapshot.size(); i++)
        if (snapshot[i].size() != snapshot[0].size()) {
            std::cerr << "Model " << i << " does not have the same size as the others. Quit!!!!\n";
            exit (0);
        }
    // Fill 1 to all the other models (models that don't have WEIGHT keyword)
    weight.resize(snapshot.size());
    for (size_t i = 0; i < weight.size(); i++)
        if (weight[i] <= 1)
            weight[i] = 1;
}

//////////////////////////////////////////////////////
// Merge H atoms into heavier atoms by distance-based
//////////////////////////////////////////////////////
void mergeH (std::vector<coordinate> &pdb_coord) {

    for (std::size_t i = 0; i < pdb_coord.size(); i++)
        if (pdb_coord[i].type == "H") {
            double min_squared = 100;
            std::size_t index;
            for (std::size_t j = 0; j < pdb_coord.size(); j++)
                if (pdb_coord[j].type != "H") {
                    double x = pdb_coord[i].x - pdb_coord[j].x;
                    double y = pdb_coord[i].y - pdb_coord[j].y;
                    double z = pdb_coord[i].z - pdb_coord[j].z;
                    double distsq = x*x + y*y + z*z;
                    if (distsq < min_squared) {
                        min_squared = distsq;
                        index = j;
                    }
                }
            pdb_coord[index].nHyd++;
        }
}

///////////////////////////////////////
// Find min and max of a pdb structure
///////////////////////////////////////
void maxmin_coord_pdb (const std::vector<coordinate> &pdb_coord,
                       coordinate &max,
                       coordinate &min) {

    double min_x = INF;    double max_x = -INF;
    double min_y = INF;    double max_y = -INF;
    double min_z = INF;    double max_z = -INF;
    for (size_t i = 0; i < pdb_coord.size(); i++) {
        if (pdb_coord[i].x < min_x)
            min_x = pdb_coord[i].x;
        else if (pdb_coord[i].x > max_x)
            max_x = pdb_coord[i].x;
        if (pdb_coord[i].y < min_y)
            min_y = pdb_coord[i].y;
        else if (pdb_coord[i].y > max_y)
            max_y = pdb_coord[i].y;
        if (pdb_coord[i].z < min_z)
            min_z = pdb_coord[i].z;
        else if (pdb_coord[i].z > max_z)
            max_z = pdb_coord[i].z;
    }
    max.x = max_x;  max.y = max_y;  max.z = max_z;
    min.x = min_x;  min.y = min_y;  min.z = min_z;
}

////////////////////////////////////////////////
// Find min and max coordinates of a trajectory
////////////////////////////////////////////////
void maxmin_coord_traj (const std::vector< std::vector<coordinate> > &snapshot,
                        coordinate &max,
                        coordinate &min) {
    coordinate max_box, min_box;
    max_box.x = -INF;   max_box.y = -INF;   max_box.z = -INF;
    min_box.x =  INF;   min_box.y =  INF;   min_box.z =  INF;

    for (size_t i = 0; i < snapshot.size(); i++) {
        coordinate max_pdb, min_pdb;
        maxmin_coord_pdb (snapshot[i], max_pdb, min_pdb);
        if (max_pdb.x > max_box.x)
            max_box.x = max_pdb.x;
        if (max_pdb.y > max_box.y)
            max_box.y = max_pdb.y;
        if (max_pdb.z > max_box.z)
            max_box.z = max_pdb.z;

        if (min_pdb.x < min_box.x)
            min_box.x = min_pdb.x;
        if (min_pdb.y < min_box.y)
            min_box.y = min_pdb.y;
        if (min_pdb.z < min_box.z)
            min_box.z = min_pdb.z;
    }
    max = max_box;    min = min_box;
}

///////////////////////////////////////////////////
// Check if this atom belongs to solute or solvent
///////////////////////////////////////////////////
bool check_solute (const coordinate &coord) {
    if ((coord.type == "Ow") or (coord.type == "Hw") or \
        (coord.type == "Na+") or (coord.type == "K+") or (coord.type == "Rb+") or (coord.type == "Cs+") or \
        (coord.type == "Mg2+") or (coord.type == "Ca2+") or (coord.type == "Sr2+") or (coord.type == "Ba2+") or \
        (coord.type == "Cl-") or (coord.type == "Br-") or (coord.type == "F-") or (coord.type == "I-"))
        return 0;
    else return 1;
}

//////////////////////////
// Get solute coordinates
//////////////////////////
std::vector<coordinate> solute_coord (const std::vector<coordinate> &pdb_coord) {
    std::vector<coordinate> solu;
    for (size_t i = 1; i < pdb_coord.size(); i++)
        if (check_solute (pdb_coord[i]))
            solu.push_back (pdb_coord[i]);
    return solu;
}

//////////////////////////////////////////////////////////////////////////
// Stripping atoms that are outside of the buffer region from the solute
//////////////////////////////////////////////////////////////////////////
std::vector<coordinate> strip_buffer_atom (const std::vector<coordinate> &model,
                                           double dcutoff,
                                           const coordinate &max_solu,
                                           const coordinate &min_solu) {
    std::vector<coordinate> solu_strip;
    for (size_t i = 0; i < model.size(); i++) {
        if (model[i].x < min_solu.x - dcutoff) continue;
        if (model[i].y < min_solu.y - dcutoff) continue;
        if (model[i].z < min_solu.z - dcutoff) continue;
        if (model[i].x > max_solu.x + dcutoff) continue;
        if (model[i].y > max_solu.y + dcutoff) continue;
        if (model[i].z > max_solu.z + dcutoff) continue;
        solu_strip.push_back (model[i]);
    }
    return solu_strip;
}

//////////////////////////////////////////////////////////////////////
// Check if an atom lies within a region or not (for smearing purpose)
///////////////////////////////////////////////////////////////////////
bool isSmear (const coordinate &atom,
              const coordinate &max_box,
              const coordinate &min_box) {
    if ((atom.x < min_box.x) or (atom.x > max_box.x)) return 0;
    else if ((atom.y < min_box.y) or (atom.y > max_box.y)) return 0;
    else if ((atom.z < min_box.z) or (atom.z > max_box.z)) return 0;
    else return 1;
}
