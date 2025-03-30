#include "opendx.hpp"
using namespace std;

namespace rism3d_c {
    opendx :: opendx(){
        // cout << "opendx object created" << endl;
    }

    opendx :: ~opendx(){
        // cout << "opendx object destroyed" << endl;
    }

    void opendx :: write_dx(string *filename, int nx, int ny, int nz, 
                            GPUtype hx, GPUtype hy, GPUtype hz, GPUtype translation[3],
                            GPUtype *u){
        std::ofstream outFile(*filename);
    
        if (!outFile) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        // Write object 1
        outFile << "object 1 class gridpositions counts " << nx << " " << ny << " " << nz << "\n";
        outFile << "origin " << -translation[0] << " " << -translation[1] << " " << -translation[2] << "\n";
        outFile << "delta " << hx << " 0.0 0.0\n";
        outFile << "delta 0.0 " << hy << " 0.0\n";
        outFile << "delta 0.0 0.0 " << hz << "\n";

        // Write object 2
        outFile << "object 2 class gridconnections counts " << nx << " " << ny << " " << nz << "\n";

        // Write object 3
        outFile << "object 3 class array type double rank 0 items " << (nx * ny * nz) << " data follows\n";

        // Write u values with fixed spacing
        const int valuesPerLine = 3;
        outFile << std::fixed << std::setprecision(6);
        int count = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    outFile << std::setw(10) << u[i * ny * (nz + 2) + j * (nz + 2) + k];
                    count++;
                    if (count % valuesPerLine == 0) {
                        outFile << "\n";
                    } else {
                        outFile << " ";
                    }
                }
            }
        }
        if (count % valuesPerLine != 0) outFile << "\n";  // Ensure final newline

        // Final object
        outFile << "object \"regular positions regular connections\" class field\n";

        outFile.close();
    }

    void opendx :: write_swapped_dx(string *filename, int nx, int ny, int nz, 
                            GPUtype hx, GPUtype hy, GPUtype hz, GPUtype translation[3],
                            GPUtype *u){
        std::ofstream outFile(*filename);

        if (!outFile) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }
    
        // Swap x and z dimensions
        outFile << "object 1 class gridpositions counts " << nz << " " << ny << " " << nx << "\n";
        outFile << "origin " << -translation[2] << " " << -translation[1] << " " << -translation[0] << "\n";
        outFile << "delta " << hz << " 0.0 0.0\n";
        outFile << "delta 0.0 " << hy << " 0.0\n";
        outFile << "delta 0.0 0.0 " << hx << "\n";
    
        // Write object 2
        outFile << "object 2 class gridconnections counts " << nz << " " << ny << " " << nx << "\n";
    
        // Write object 3
        outFile << "object 3 class array type double rank 0 items " << (nx * ny * nz) << " data follows\n";
    
        // Write u values with fixed spacing
        const int valuesPerLine = 3;
        outFile << std::fixed << std::setprecision(6);
        int count = 0;
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    outFile << std::setw(10) << u[i * ny * (nz + 2) + j * (nz + 2) + k];
                    count++;
                    if (count % valuesPerLine == 0) {
                        outFile << "\n";
                    } else {
                        outFile << " ";
                    }
                }
            }
        }
        if (count % valuesPerLine != 0) outFile << "\n";  // Ensure final newline
    
        // Final object
        outFile << "object \"regular positions regular connections\" class field\n";
    
        outFile.close();
    }
}