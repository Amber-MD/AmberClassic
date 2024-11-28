/*

The MoFT software found here is copyright (c) 2019 by
George M. Giambasu, Darrin M. York and David A. Case.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License in the
LICENSE file.  If not, see <http://www.gnu.org/licenses/>.

Users of the MoFT suite of programs found here are requested to acknowledge
use of the software in reports and publications.  Such acknowledgement
should include the following citations:

(1) "Ion counting from explicit-solvent simulations and 3D-RISM"
GM Giambaşu, T Luchko, D Herschlag, DM York, DA Case
Biophysical Journal 106 (4), 883-894 doi:10.1016/j.bpj.2014.01.021

(2) "Competitive interaction of monovalent cations with DNA from 3D-RISM"
GM Giambaşu, MK Gebala, MT Panteva, T Luchko, DA Case, DM York
Nucleic acids research 43 (17), 8405-8415 doi:10.1093/nar/gkv830

(3) "Predicting site-binding modes of ions and water to nucleic acids using molecular solvation theory"
GM Giambasu, DA Case, DM York
Journal of the American Chemical Society doi:10.1021/jacs.8b11474

*/

#ifndef DENSITY_HPP_
#define DENSITY_HPP_

#include <algorithm>
#include <array>
#include "boosthead.hpp"
#include "readdx.hpp"
#include "readccp4.hpp"
#include "scattering.hpp"
#include "utilities/utilities.hpp"
#include "metatwist/lutilities.hpp"
#include "metaFFT/fftwwrap.hpp"
#include "metaFFT/ConvCorr.hpp"



template <typename T>
class Density{
    
public:
    // List of defined types:
    typedef boost::multi_array<T, 3> Tarray;
    typedef ublas::matrix<T> Tmat;
    typedef ublas::vector<T> Tvec;
    
    // Constructors
    
    // Constructs a density reading from list of files.
    Density(const std::vector <std::string> & filenames, const std::vector<std::string> & cspecies );
    
    // Get the dataset name
    std::string getName(){return this->dataname;}

    // Set data of grid point.
    void set(size_t i, size_t j, size_t k,const T & data){ this->data3d[i][j][k] = data;}
    
    // Set the entire grid.
    void set(const Tarray & data);
    
    // Set, get rotation matrix.
    void setrotm(const Tmat & atmat);
    Tmat getrotm(){return this->rotm;}
    
    // Set, get origin.
    void setorig(const Tvec & avec);
    Tvec getorig(){return this->orig;}
    
    
    // set, get chemical species
    
    void setspecies(const std::string & spc){this->species = spc;}

    void setcharge(const T& charge){this->partialcharge=charge;}
    std::string getspecies(){return this->species;}
    
    
    // Get grid size.
    ublas::vector<size_t> getGridSize(){
        ublas::vector<size_t> dsizes(3,0);
        dsizes(0) = this->data3d.shape()[0]; dsizes(1) = this->data3d.shape()[1]; dsizes(2) = this->data3d.shape()[2];
        return dsizes;
    }
    
    void setBulkC(const T & b){this->bulkc = b * 6.023 * 1e-4;}

    void average (){
         // taking an average in the case of multiple input densities
         if (count > 1) {
             for (size_t i = 0; i != this->data3d.num_elements(); ++i) {
                 this->data3d.data()[i] = this->data3d.data()[i]/count;
             }
         }
    }
    
    // Gets density value given x,y,z using tri-linear interpolation.
    // assumes 0 outside the box
    T get(const Tvec & xyz);
    // assumes a periodic box
    T getper(const Tvec & xyz);

    // Returns a 3d-vector with the density grid origin.
    Tvec getOrigin(){ return this->orig;}
    
    // Computes excess number of particles.
    T getIntExcess(T & bulkdens);
    
    // Say whether the density object has any actual 3d-data read in.
    bool getreadany(){ return this->readany;}
    
    // Project the entire grid on X, Y, Z axes
    void getProjXYZ(Tvec & projx, Tvec & projy, Tvec & projz, const T & bulkdens);
    
    // Convolute data using different kernels
    void convolute(const T & sigma, const size_t convtype = 3);
    void convolutex(const T & sigma, const size_t convtype = 3);
    
    // Compute Laplacian.
    void getLaplacian(Density<T> & LDensity);
    
    // Compute Normalize Laplacian.
    void getNLaplacian(Density<T> & LDensity);
    
    void getNpLaplacian(Density<T> & LDensity);
    
    // take the negative log of the data, overwrites initial data.
    void nlogtrans();
    
    // get electron density based on atomic structure factors
    
    void getRhoEl();
    void getRhoElReal();


    // get charge density

    void getQDensity(Density<T> & QDensity);
    // cut the reciprocal space resolutions

    void cutResolution(const T & min, const T & max);
    
    // Carry out blob analysis based on the sign of the Laplacian.
    void blobs(Density<T> & laplacian, const T & threshold, const T & bulkc = 1.00);
    void blobs_periodic(Density<T> & laplacian, const T & threshold);
    
    // Write the density to a dx, ccp4, or mrc file.
    void writedxfile(const std::string &filenameout);
    void writeccp4(const std::string &filenameout);

    // print details
    void printDetails();
    
private:
    // array of data
    Tarray data3d;
    
    // rotation matrix
    Tmat rotm,irotm;
    
    // origin
    Tvec orig;
    
    // whether or not any density has been read
    bool readany;
    
    // name of the dataset
    
    std::string dataname, species, residue;
    
    std::string chemicalspecies;
    T partialcharge;
    
    T Vvox, Vol, bulkc;

    // real and reciprocal vectors
    Tvec recipra,reciprb,reciprc;
    Tvec a,b,c;

    // metric tensors

    Tmat metrictensor, imetrictensor;

    size_t count;
    
    
    
};


template <typename T>
Density<T>::Density(const std::vector <std::string> & filenames , const std::vector<std::string> & cspecies) {
    
    this->readany = false;
    this->count = 0;
    this->partialcharge = 1.0;


    for (const std::string & filename : filenames){
        
        std::vector<std::string> stokens;
        stokens =  boost::copy_range<std::vector<std::string>>(filename |boost::adaptors::tokenized( std::string("\\."), -1 ) );
        std::string ext,ext1;
        ext = stokens[stokens.size()-1];
        ext1 = stokens[stokens.size()-2];
        

        if( ext == "dx") {
            if (util::files::readdx(filename, this->orig,this->rotm,this->data3d)){
                //std::cout << "# >> Read density from: " << filename << std::endl ;
                this->readany = true;
                count++;
            }
        } else if ( ext=="bz2")
        {
            
            if (ext1=="dx"){
                
                if (util::files::readdxbz2(filename, this->orig,this->rotm,this->data3d)){
                    //std::cout << "# >> Read density from: " << filename << std::endl ;
                    this->readany = true;
                    count++;
                }
            }
            
        }
        else if ( ext=="gz")
        {
            if (ext1=="dx") {
                if (util::files::readdxgz(filename, this->orig,this->rotm,this->data3d)){
                    // std::cout << "# >> Read density from: " << filename << std::endl ;
                    this->readany = true;
                    count++;
                }
            }
            
            if (ext1=="xyzv"){
                
                if (util::files::readxyzvgz(filename, this->orig,this->rotm,this->data3d)){
                    //std::cout << "# >> Read density from: " << filename << std::endl ;
                    this->readany = true;
                    count++;
                }
                
            }
        }
        
        else if ( ext=="xyzv")
        {
            
            if (util::files::readxyzv(filename, this->orig,this->rotm,this->data3d)){
                //std::cout << "# >> Read density from: " << filename << std::endl ;
                this->readany = true;
                count++;
            }
            
        }
        
        else if ( ext=="ccp4" || ext=="mrc" )
        {

            if (util::files::readccp4(filename, this->orig,this->rotm,this->data3d)){
                //std::cout << "# >> Read density from: " << filename << std::endl ;
                this->readany = true;
                count++;
            }
        }
        
        else {
            std::cout << "Wrong extension!" << std::endl;
            
        }
        
    }
    
    
    
    if (this->readany == false){
        std::cout << "# Density data is empty!" << std::endl;
    }
    
    
    
 //   // taking an average in the case of multiple input densities
 //   if (count > 1) {
 //
 //       // std::cout << "# Averaging out the input densities!" << std::endl;
 //
 //       for (size_t i = 0; i != this->data3d.num_elements(); ++i) {
 //
 //           this->data3d.data()[i] = this->data3d.data()[i]/count;
 //       }
 //   }
    
    if (this->readany == true) {


        if(!(util::InvertMatrix(rotm,irotm))){std::cout << "#>> Unable to inverted rot matrix." << std::endl ;};


        // compute voxel volume
        
        typedef typename Tarray::index index;
        index xdim = this->data3d.shape()[0];
        index ydim = this->data3d.shape()[1];
        index zdim = this->data3d.shape()[2];
        
        this->species = cspecies[0];
        
        if (cspecies.size()>1) {
            this->residue = cspecies[1];
        } else {
            this->residue = this->species;
        }
        
        
        this->bulkc = 1.00;
        
        this->Vvox = std::abs((-this->rotm(0,2)*this->rotm(1,1) + this->rotm(0,1)*this->rotm(1,2))*this->rotm(2,0) +
                              ( this->rotm(0,2)*this->rotm(1,0) - this->rotm(0,0)*this->rotm(1,2))*this->rotm(2,1) +
                              (-this->rotm(0,1)*this->rotm(1,0) + this->rotm(0,0)*this->rotm(1,1))*this->rotm(2,2))   ;
        
        // cell volume
        this->Vol = this->Vvox*xdim*ydim*zdim;
        
        // compute reciprocal vectors
        
        this->recipra.resize(3);
        this->reciprb.resize(3);
        this->reciprc.resize(3);
        this->recipra(0)=(-this->rotm(1,2)*this->rotm(2,1) + this->rotm(1,1)*this->rotm(2,2))*ydim*zdim/this->Vol;
        this->recipra(1)=( this->rotm(1,2)*this->rotm(2,0) - this->rotm(1,0)*this->rotm(2,2))*ydim*zdim/this->Vol;
        this->recipra(2)=(-this->rotm(1,1)*this->rotm(2,0) + this->rotm(1,0)*this->rotm(2,1))*ydim*zdim/this->Vol;
        
        this->reciprb(0)=( this->rotm(0,2)*this->rotm(2,1) - this->rotm(0,1)*this->rotm(2,2))*xdim*zdim/this->Vol;
        this->reciprb(1)=(-this->rotm(0,2)*this->rotm(2,0) + this->rotm(0,0)*this->rotm(2,2))*xdim*zdim/this->Vol;
        this->reciprb(2)=( this->rotm(0,1)*this->rotm(2,0) - this->rotm(0,0)*this->rotm(2,1))*xdim*zdim/this->Vol;
        
        this->reciprc(0)=(-this->rotm(0,2)*this->rotm(1,1) + this->rotm(0,1)*this->rotm(1,2))*xdim*ydim/this->Vol;
        this->reciprc(1)=( this->rotm(0,2)*this->rotm(1,0) - this->rotm(0,0)*this->rotm(1,2))*xdim*ydim/this->Vol;
        this->reciprc(2)=(-this->rotm(0,1)*this->rotm(1,0) + this->rotm(0,0)*this->rotm(1,1))*xdim*ydim/this->Vol;




        // compute metric tensors
        this->metrictensor.resize(3,3);
        this->imetrictensor.resize(3,3);

       // this->imetrictensor * = 0;


        //Tmat imetrictensor(3,3,0);
        this->imetrictensor = this->rotm;
        Tmat im (3,3,0); im(0,0) = xdim; im(1,1) = ydim; im(2,2) = zdim;
        im = ublas::prod(imetrictensor,im);
        im = ublas::prod(im, ublas::trans(im));
        this->metrictensor = im;
        util::InvertMatrix(im,this->imetrictensor);


        
        
        // compute real vectors

        Tmat tmp (3,3,0);
        tmp(0,0) = xdim; tmp(1,1) = ydim; tmp(2,2) = zdim;
        tmp =  ublas::prod(this->rotm,tmp);
        this->a.resize(3); this->a(0) = tmp(0,0); this->a(1) = tmp(1,0); this->a(2) = tmp(2,0);
        this->b.resize(3); this->b(0) = tmp(0,1); this->b(1) = tmp(1,1); this->b(2) = tmp(2,1);
        this->c.resize(3); this->c(0) = tmp(0,2); this->c(1) = tmp(1,2); this->c(2) = tmp(2,2);

        // a, b, c, angles
        
        
        
        std::vector<std::string> stokens;
        stokens =  boost::copy_range<std::vector<std::string>>
        (filenames[0] |boost::adaptors::tokenized( std::string("(\\.dx)|(\\.ccp4)|(\\.mrc)|(\\.xyzv)"), -1 ) );
        this->dataname = stokens[0];
        std::cout << boost::format("#\n#\n# Volumetric data initiated.\n# name        : %s" ) % this->dataname << std::endl;
        for (size_t i = 0 ; i!= filenames.size(); ++i){
            std::cout << boost::format("# file        : %s") % filenames[i]  << std::endl;
        }
        std::cout << "#\n#\n#" << std::endl;
    }
}



template <typename T>
void Density<T>::printDetails(){
    std::cout << boost::format("#\n#\n# Volumetric data details.\n# name        : %s" ) % this->dataname << std::endl;
    std::cout << boost::format("# species     : %s ") % this->species << std::endl;
    std::cout << boost::format("# residue     : %s ") % this->residue << std::endl;
    std::cout << boost::format("# bulk dens(particle * A^-3): %8.3E ") % this->bulkc << std::endl;
    std::cout << boost::format("# bulk conc(M): %8.3f ") % (this->bulkc/(6.023*1e-4)) << std::endl;
    std::cout << "# " << std::endl;
    std::cout << boost::format("# origin [A]  : %8.3f %8.3f %8.3f") % this->orig(0) % this->orig(1) % this->orig(2) << std::endl;
    std::cout << boost::format("# skew        : %8.3f %8.3f %8.3f") % this->rotm(0,0) % this->rotm(0,1) % this->rotm(0,2) << std::endl;
    std::cout << boost::format("# skew        : %8.3f %8.3f %8.3f") % this->rotm(1,0) % this->rotm(1,1) % this->rotm(1,2) << std::endl;
    std::cout << boost::format("# skew        : %8.3f %8.3f %8.3f") % this->rotm(2,0) % this->rotm(2,1) % this->rotm(2,2) << std::endl;
    std::cout << boost::format("# max value   : %8.3f ") % *std::max_element(this->data3d.data(),this->data3d.data()+this->data3d.num_elements()) << std::endl;
    std::cout << boost::format("# min value   : %8.3f ") % *std::min_element(this->data3d.data(),this->data3d.data()+this->data3d.num_elements()) << std::endl;
    std::cout << boost::format("# grid size   : %i x %i x %i ") % this->data3d.shape()[0] % this->data3d.shape()[1] % this->data3d.shape()[2] << std::endl;
    std::cout << "# " << std::endl;
    std::cout << boost::format("# a  [A]      : %8.4f %8.4f %8.4f") % this->a(0) % this->a(1) % this->a(2) << std::endl;
    std::cout << boost::format("# b  [A]      : %8.4f %8.4f %8.4f") % this->b(0) % this->b(1) % this->b(2) << std::endl;
    std::cout << boost::format("# c  [A]      : %8.4f %8.4f %8.4f") % this->c(0) % this->c(1) % this->c(2) << std::endl;
    std::cout << boost::format("# a* [A^-1]   : %8.4f %8.4f %8.4f") % this->recipra(0) % this->recipra(1) % this->recipra(2) << std::endl;
    std::cout << boost::format("# b* [A^-1]   : %8.4f %8.4f %8.4f") % this->reciprb(0) % this->reciprb(1) % this->reciprb(2) << std::endl;
    std::cout << boost::format("# c* [A^-1]   : %8.4f %8.4f %8.4f") % this->reciprc(0) % this->reciprc(1) % this->reciprc(2) << std::endl;
    std::cout << "# " << std::endl;
    std::cout << boost::format("# Volume [A^3]: %12.3f ") % this->Vol << std::endl;
    std::cout << boost::format("# Vvoxel [A^3]: %12.3f ") % this->Vvox << std::endl;
}


template <typename T>
T Density<T>::getper(const Tvec & xyz){

  Tvec v = ublas::prod(this->irotm,xyz - this->orig);
  int v1min = std::floor(v(0));
  int v2min = std::floor(v(1));
  int v3min = std::floor(v(2));
  int v1max = v1min+1;
  int v2max = v2min+1;
  int v3max = v3min+1;

  int s0 = this->data3d.shape()[0];
  int s1 = this->data3d.shape()[1];
  int s2 = this->data3d.shape()[2];

//  v1min = (v1min>=0)? v1min % s0 : ( s0 + v1min % s0 ) % s0;
//  v2min = (v2min>=0)? v2min % s1 : ( s1 + v2min % s1 ) % s1;
//  v3min = (v3min>=0)? v3min % s2 : ( s2 + v3min % s2 ) % s2;
//
//  v1max = (v1max>=0)? v1max % s0 : ( s0 + v1max % s0 ) % s0;
//  v2max = (v2max>=0)? v2max % s1 : ( s1 + v2max % s1 ) % s1;
//  v3max = (v3max>=0)? v3max % s2 : ( s2 + v3max % s2 ) % s2;


  v1min = (v1min>=0)? v1min % s0 : (5*s0 + v1min) % s0;
  v2min = (v2min>=0)? v2min % s1 : (5*s1 + v2min) % s1;
  v3min = (v3min>=0)? v3min % s2 : (5*s2 + v3min) % s2;

  v1max = (v1max>=0)? v1max % s0 : (5*s0 + v1max) % s0;
  v2max = (v2max>=0)? v2max % s1 : (5*s1 + v2max) % s1;
  v3max = (v3max>=0)? v3max % s2 : (5*s2 + v3max) % s2;


 //   std::cout << " v= " <<   v1min ;
 //   std::cout << " v= " <<   v2min ;
 //   std::cout << " v= " <<   v3min ;
 //   std::cout << " v= " <<   v1max ;
 //   std::cout << " v= " <<   v2max ;
 //   std::cout << " v= " <<   v3max << std::endl;

    return util::blend103( v(0) - T(std::floor(v(0))),
                           v(1) - T(std::floor(v(1))),
                           v(2) - T(std::floor(v(2))),
                           this->data3d[v1min][v2min][v3min],
                           this->data3d[v1min][v2min][v3max],
                           this->data3d[v1min][v2max][v3min],
                           this->data3d[v1min][v2max][v3max],
                           this->data3d[v1max][v2min][v3min],
                           this->data3d[v1max][v2min][v3max],
                           this->data3d[v1max][v2max][v3min],
                           this->data3d[v1max][v2max][v3max]
                         );

}







template <typename T>
T Density<T>::get(const Tvec & xyz){
    
    //Tvec v = xyz - this->orig;
    Tvec v = ublas::prod(this->irotm,xyz - this->orig);
    
    if ( v(0)>=0.0 && v(1)>=0.0 && v(2)>=0.0) {
        
        size_t v1min = std::floor(v(0));
        size_t v2min = std::floor(v(1));
        size_t v3min = std::floor(v(2));
        size_t v1max = v1min+1;
        size_t v2max = v2min+1;
        size_t v3max = v3min+1;

        if (v1max < this->data3d.shape()[0]
            && v2max < this->data3d.shape()[1]
            && v3max < this->data3d.shape()[2]
            )
        {
            T r = v(0) - v1min;
            T s = v(1) - v2min;
            T t = v(2) - v3min;
            T x000 = this->data3d[v1min][v2min][v3min];
            T x100 = this->data3d[v1max][v2min][v3min];
            T x110 = this->data3d[v1max][v2max][v3min];
            T x010 = this->data3d[v1min][v2max][v3min];
            T x001 = this->data3d[v1min][v2min][v3max];
            T x101 = this->data3d[v1max][v2min][v3max];
            T x111 = this->data3d[v1max][v2max][v3max];
            T x011 = this->data3d[v1min][v2max][v3max];
            
            return util::blend103(r,s,t,x000, x001, x010, x011, x100, x101, x110, x111);
        } else {
            //std::cout <<"#>> Out of bounds!" << xyz << std::endl;
            return 0.0;
        }
        
    }
    else {
        
        
        return 0;
    }
}

template <typename T>
T Density<T>::getIntExcess(T & bulkdens){
    T intgrl = 0.0;
    std::cout << "# >> Bulk density: " << bulkdens << std::endl;
    for (size_t i = 0 ; i != data3d.num_elements() ; ++i ) {
        intgrl += (this->data3d.data()[i]-bulkdens);
    }
    return (intgrl * this->Vvox);
}

template <typename T>
void Density<T>::getProjXYZ(Tvec & projx, Tvec & projy, Tvec & projz, const T & bulkdens){
    
    // resize projx,y,x and set to zero
    size_t xdim = this->data3d.shape()[0];
    size_t ydim = this->data3d.shape()[1];
    size_t zdim = this->data3d.shape()[2];
    
    projx.resize(xdim);
    projy.resize(ydim);
    projz.resize(zdim);
    
    
    for (size_t i = 0; i != projx.size();++i) { projx(i) =0.0;};
    for (size_t i = 0; i != projy.size();++i) { projy(i) =0.0;};
    for (size_t i = 0; i != projz.size();++i) { projz(i) =0.0;};
    T sum = 0.0;
    
    
    
    // loop over all points, accumulate data, average
    for (size_t i = 0; i != xdim; ++i){
        for (size_t j = 0; j != ydim; ++j){
            for (size_t k = 0; k != zdim; ++k){
                T tmp = (this->data3d[i][j][k] - bulkdens);
                tmp *= this->Vvox;
                projz(k) += tmp;
                projy(j) += tmp;
                projx(i) += tmp;
                sum += tmp;
                
            }
        }
    }
    
    std::cout << "# Sum  : " <<  std::setprecision (15) <<  sum << std::endl;
    std::cout << "# Sum x: " <<  std::setprecision (15) <<  ublas::sum(projx) << " " << projx.size() << std::endl;
    std::cout << "# Sum y: " <<  std::setprecision (15) <<  ublas::sum(projy) << " " << projy.size() << std::endl;
    std::cout << "# Sum z: " <<  std::setprecision (15) <<  ublas::sum(projz) << " " << projz.size() << std::endl;
    
    
    
    
}


template <typename T>
void Density<T>::convolutex(const T & sigma, const size_t convtype){
    
    typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
    //typedef boost::multi_array<double, 3> TArray;
    typedef typename Tarray::index index;
    index xdim = this->data3d.shape()[0];
    index ydim = this->data3d.shape()[1];
    index zdim = this->data3d.shape()[2];
    
    boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
    boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, (zdim/2)+1}};
    index dh = (zdim/2)-1;

    index dhx,dhy,dhz;


    if (xdim % 2 == 0) dhx = (xdim/2)-1; else { dhx = ((xdim-1)/2) -1;}
    if (ydim % 2 == 0) dhy = (ydim/2)-1; else { dhy = ((ydim-1)/2) -1;}
    if (zdim % 2 == 0) dhz = (zdim/2)-1; else { dhz = ((zdim-1)/2) -1;}

    ublas::vector<T> centr (3,0);
    centr(0)=dhx; centr(1)=dhy; centr(2)=dhz;
    //std::cout << " central location" << centr << std::endl;
    centr = ublas::prod(this->rotm,centr);
    //std::cout << " central location" << centr << std::endl;
    
    // where complex fft are stored
    Tarray_c  Afft(shapeft), Bfft(shapeft);
    // where real data is stored
    Tarray A(shape), B(shape);
    
    
    // init to zero
    for (int i = 0 ; i != xdim; ++i){
        for (int j = 0; j!=ydim; ++j){
            for (int k = 0; k != zdim; ++k){
                A[i][j][k] = 0.0;
                B[i][j][k] = 0.0;
            }
        }
    }
    
    T sum = 0.0;
    
    if (convtype == 1){
        std::cout << "# Gaussian Filter." << std::endl;
        
        //double sum = 0.0;
        double tsq=2*sigma*sigma;
        for (int i = 0 ; i != xdim; ++i){
            for (int j = 0; j!=ydim; ++j){
                for (int k = 0; k != zdim; ++k){
                    Tvec point (3,0);
                    point (0) = i; point (1) = j; point (2) = k;
                    point = ublas::prod(this->rotm,point);
                    point = point - centr;
                    B[i][j][k] = std::exp(-1*ublas::inner_prod(point,point)/tsq);
                    sum += B[i][j][k];
                }
            }
        }
        
        std::cout << sum << std::endl;
        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
        
    } else if (convtype ==2) {
        
        std::cout << "# Box Filter." << std::endl;
        
        //        int lmin = dh - (int(sigma)), lmax = dh+(int(sigma)) +1;
        
        //        for (int i = lmin ; i != lmax; ++i){
        //            for (int j = lmin; j!= lmax; ++j){
        //                for (int k = lmin; k != lmax; ++k){
        //                    B[i][j][k] = 1.0  ;
        //                    sum += B[i][j][k];
        //                }
        //            }
        //        }
        //
        
        //        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
        
        //    }
        //    else if (convtype==3) {
        //        std::cout << "# Sinc Filter." << std::endl;
        //        // hardcode the order
        //        int lsc = 3;
        //
        //
        //        int a = int(sigma+1);
        //
        //        int lmin = dh-a*lsc, lmax = dh+a*lsc+1 ;
        //
        //
        //
        //        for (int i = lmin ; i != lmax; ++i){
        //            for (int j = lmin; j!= lmax; ++j){
        //                for (int k = lmin; k != lmax; ++k){
        //
        //                    double xd =std::fabs((i-dh)*M_PI/a);
        //                    double yd =std::fabs((j-dh)*M_PI/a);
        //                    double zd =std::fabs((k-dh)*M_PI/a);
        //
        //                    double tmp = boost::math::sinc_pi(xd) * boost::math::sinc_pi(xd/lsc)
        //                    * boost::math::sinc_pi(yd) * boost::math::sinc_pi(yd/lsc)
        //                    * boost::math::sinc_pi(zd) * boost::math::sinc_pi(zd/lsc);
        //
        //                    B[i][j][k] = tmp;
        //                    sum += tmp;
        //
        //                }
        //            }
        //        }
        //        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
    } else if (convtype==4) {
        
        std::cout << "# Laplacian of a Gaussian Filter." << std::endl;
        double tsq=(2*sigma*sigma);
        double scale = 1.0/(2.0 * std::sqrt(2.0) * std::pow(M_PI,1.5) * std::pow(sigma,7));
        for (int i = 0 ; i != xdim; ++i){
            for (int j = 0; j!=ydim; ++j){
                for (int k = 0; k != zdim; ++k){
                    Tvec point (3,0);
                    point (0) = i; point (1) = j; point (2) = k;
                    point = ublas::prod(this->rotm,point);
                    point = point - centr;
                    T r = ublas::inner_prod(point, point);
                    B[i][j][k] = std::exp(-r/tsq) * ( -1.5*tsq + r)*scale ;
                }
            }
        }
        
        
    } else {
        
        std::cout << "# This filter has not been implemented yet." << std::endl;
        
        //return 1;
    }
    
    
    
    
    
    

    
    // A is the delta function, used to shift filter to zero
    
    //A[xdim-dh][ydim-dh][zdim-dh] = 1.0;
    A[xdim-dhx][ydim-dhy][zdim-dhz] = 1.0;
    // inplace convolution, overwrite B
    metafft::convolution(B,A);
    
    // inplace convolution, overwrite data3d
    metafft::convolution(this->data3d,B);
    
   // for (size_t i = 0; i!= data3d.num_elements();++i){
   //
   //     this->data3d.data()[i] *= this->bulkc;
   // }
    
    
    
    
    
    
    
}


/*
template <typename T>
void Density<T>::convolute(const T & sigma, const size_t convtype){
    
    typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
    //typedef boost::multi_array<double, 3> TArray;
    typedef typename Tarray::index index;
    index xdim = this->data3d.shape()[0];
    index ydim = this->data3d.shape()[1];
    index zdim = this->data3d.shape()[2];
    
    boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
    boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, (zdim/2)+1}};
    index dh = (zdim/2)-1;
    
    // where complex fft are stored
    Tarray_c  Afft(shapeft), Bfft(shapeft);
    // where real data is stored
    Tarray A(shape), B(shape) ;
    
    // init to zero
    for (int i = 0 ; i != xdim; ++i){
        for (int j = 0; j!=ydim; ++j){
            for (int k = 0; k != zdim; ++k){
                A[i][j][k] = 0.0;
                B[i][j][k] = 0.0;
            }
        }
    }
    
    double sum = 0.0;
    
    if (convtype == 1){
        std::cout << "# Gaussian Filter." << std::endl;
        
        double sum = 0.0;
        double tsq=2*sigma*sigma;
        for (int i = 0 ; i != xdim; ++i){
            for (int j = 0; j!=ydim; ++j){
                for (int k = 0; k != zdim; ++k){
                    double xd = i-dh;
                    double yd = j-dh;
                    double zd = k-dh;
                    
                    
                    
                    
                    
                    B[i][j][k] = std::exp(-xd*xd/tsq) * std::exp(-yd*yd/tsq) * std::exp(-zd*zd/tsq)  ;
                    sum += B[i][j][k];
                }
            }
        }
        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
        
        
    } else if (convtype ==2) {
        
        std::cout << "# Box Filter." << std::endl;
        
        int lmin = dh - (int(sigma)), lmax = dh+(int(sigma)) +1;
        
        for (int i = lmin ; i != lmax; ++i){
            for (int j = lmin; j!= lmax; ++j){
                for (int k = lmin; k != lmax; ++k){
                    B[i][j][k] = 1.0  ;
                    sum += B[i][j][k];
                }
            }
        }
        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
        
    } else if (convtype==3) {
        std::cout << "# Sinc Filter." << std::endl;
        // hardcode the order
        int lsc = 3;
        
        
        int a = int(sigma+1);
        
        int lmin = dh-a*lsc, lmax = dh+a*lsc+1 ;
        
        
        
        for (int i = lmin ; i != lmax; ++i){
            for (int j = lmin; j!= lmax; ++j){
                for (int k = lmin; k != lmax; ++k){
                    
                    double xd =std::fabs((i-dh)*M_PI/a);
                    double yd =std::fabs((j-dh)*M_PI/a);
                    double zd =std::fabs((k-dh)*M_PI/a);
                    
                    double tmp = boost::math::sinc_pi(xd) * boost::math::sinc_pi(xd/lsc)
                    * boost::math::sinc_pi(yd) * boost::math::sinc_pi(yd/lsc)
                    * boost::math::sinc_pi(zd) * boost::math::sinc_pi(zd/lsc);
                    
                    B[i][j][k] = tmp;
                    sum += tmp;
                    
                }
            }
        }
        for (size_t i = 0; i != B.num_elements();++i) { B.data()[i] /= sum; }
    } else if (convtype==4) {
        
        std::cout << "# Laplacian of a Gaussian Filter." << std::endl;
        double tsq=2*sigma*sigma;
        double scale = 2.0 * std::sqrt(2.0) * std::pow(M_PI,1.5) * std::pow(sigma,7);
        for (int i = 0 ; i != xdim; ++i){
            for (int j = 0; j!=ydim; ++j){
                for (int k = 0; k != zdim; ++k){
                    double xd = i-dh;
                    double yd = j-dh;
                    double zd = k-dh;
                    xd *=xd;
                    yd *=yd;
                    zd *=zd;
                    
                    xd += yd;
                    xd += zd;
                    
                    B[i][j][k] = std::exp(-xd/tsq) * ( -1.5*tsq + xd) /scale ;
                    //sum += B[i][j][k];
                }
            }
        }
        
        
    }
    
    
    
    // A is the delta function, used to shift filter to zero
    A[xdim-dh][ydim-dh][zdim-dh] = 1.0;
    
    // inplace convolution, overwrite B
    metafft::convolution(B,A);
    
    // inplace convolution, overwrite data3d
    metafft::convolution(this->data3d,B);
    
}
*/



template <typename T>
void Density<T>::getRhoElReal(){



  std::cout << "#\n#\n# Computing electron distribution" << std::endl;
  std::cout << "#    Particle Distribution: " << this->dataname << std::endl;


  if (scattering::themap.count(this->species) > 0 ){

      scattering::scatterer<T,12> myscatterer (scattering::themap[this->species]);
      std::cout << "#    Atom form factors for: " << this->species << " having "
                <<  int(std::round(myscatterer.getf(T(0)))) << " electrons."<<   std::endl;

      std::cout << "#    Total number of particles: "<<
          std::accumulate(this->data3d.data(), this->data3d.data()+this->data3d.num_elements(),0.0) * this->bulkc * this->Vvox
                << std::endl;

      typedef typename Tarray::index index;
      index xdim = this->data3d.shape()[0];
      index ydim = this->data3d.shape()[1];
      index zdim = this->data3d.shape()[2];

      boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};

      index dhx,dhy,dhz;


      if (xdim % 2 == 0) dhx = (xdim/2)-1; else { dhx = ((xdim-1)/2) -1;}
      if (ydim % 2 == 0) dhy = (ydim/2)-1; else { dhy = ((ydim-1)/2) -1;}
      if (zdim % 2 == 0) dhz = (zdim/2)-1; else { dhz = ((zdim-1)/2) -1;}

      ublas::vector<T> centr (3,0);
      centr(0)=dhx; centr(1)=dhy; centr(2)=dhz;

      centr = ublas::prod(this->rotm,centr);

      Tarray A(shape), B(shape);

      for (int i = 0 ; i != xdim; ++i){
          for (int j = 0; j!=ydim; ++j){
              for (int k = 0; k != zdim; ++k){
                  Tvec point (3,0);
                  point (0) = i; point (1) = j; point (2) = k;
                  point = ublas::prod(this->rotm,point);
                  point = point - centr;
                  B[i][j][k] = myscatterer.getrhor(ublas::inner_prod(point, point));
                  A[i][j][k] = 0;
                }
            }
        }

      T numel = std::accumulate(B.data(), B.data()+B.num_elements(),T(0.0)) ;
      //std::cout << "#    Total number of electrons in kernel: "<< numel * this->Vvox << std::endl;

      for (size_t i = 0; i != B.num_elements(); ++i)  B.data()[i] *= (myscatterer.getf(0.0)/numel)  ;


      numel = std::accumulate(B.data(), B.data()+B.num_elements(),T(0.0)) ;

      // A is the delta function, used to shift filter to zero
      A[xdim-dhx][ydim-dhy][zdim-dhz] = 1.0;


      // inplace convolution, overwrite B, shift to origin
      metafft::convolution(B,A);

      metafft::convolution(this->data3d,B);

      for (size_t i = 0; i != this->data3d.num_elements(); ++i)  this->data3d.data()[i] *= this->bulkc ;

      std::cout << "#    Total number of electrons: "<<
          std::accumulate(this->data3d.data(), this->data3d.data()+this->data3d.num_elements(),0.0) * this->Vvox
                << std::endl;



    } else {

      std::cout << "# The atomic scattering factors for the specified chemical species do not exist." << std::endl;
      std::cout << "# Consider choosing from this list: " << std::endl;

      for (auto const& cspec : scattering::themap) {
          std::cout << "# " << cspec.first << std::endl;
        }

    }

}

template <typename T>
void Density<T>::getRhoEl( ){


  std::cout << "#\n#\n# Computing electron density" << std::endl;
  std::cout << "#    Particle Distribution: " << this->dataname << std::endl;


  if (scattering::themap.count(this->species) > 0 ){

      scattering::scatterer<T,12> myscatterer (scattering::themap[this->species]);

      std::cout << "#    Atom form factors for: "
                << this->species << " having  a total of "
                <<  int(std::round(myscatterer.getf(T(0)))) << " electrons."
                <<   std::endl;



      typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
      //typedef boost::multi_array<double, 3> TArray;
      typedef typename Tarray::index index;
      index xdim = this->data3d.shape()[0];
      index ydim = this->data3d.shape()[1];
      index zdim = this->data3d.shape()[2];

      boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
      boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, zdim}};


      // where complex fft are stored
      Tarray_c  Afft(shapeft);
      // where real data is stored
      Tarray A(shape);

      // copy data3d into A
      T integral = 0;
      for (int i = 0 ; i != xdim; ++i){
          for (int j = 0; j!=ydim; ++j){
              for (int k = 0; k != zdim; ++k){
                  A[i][j][k] = this->data3d[i][j][k] * this->bulkc;
                  integral += A[i][j][k] * this->Vvox;
                }
            }
        }

      std::cout << "#   integral over initial particle distribution: " << integral << " particles." << std::endl;

      // get dens3d ft
      metafft::FFT(A,Afft,-1);

      for (int i = -xdim/2 ; i != xdim/2; ++i){
          int ii = ((i<0)? i+xdim : i );
          for (int j = -ydim/2; j != ydim/2; ++j){
               int jj = ((j<0)? j+ydim : j);
              for (int k = -zdim/2;  k != zdim/2; ++k){
                   int kk = ((k<0)? k+zdim : k);
                       Tvec v(3); v(0) = i; v(1) = j; v(2) = k;
                       Afft[ii][jj][kk] *= myscatterer.getfc(  ublas::inner_prod(ublas::prod(v,this->imetrictensor),v) ) ;
                }
            }
        }


      metafft::FFT(this->data3d,Afft,1);

      integral = 0;

      for (size_t i = 0 ; i!=this->data3d.num_elements(); ++i) {
          this->data3d.data()[i] /= this->data3d.num_elements();
          integral += this->data3d.data()[i];
        }

      std::cout << "#    Integral over resulting electron density: " << integral * this->Vvox << " electrons." << std::endl;

    } else {

      std::cout << "#    The atomic scattering factors for the specified chemical species do not exist." << std::endl;
      std::cout << "#    Consider choosing from this list: " << std::endl;

      for (auto const& cspec : scattering::themap) {
          std::cout << "# " << cspec.first << std::endl;
        }

    }
    
}


template <typename T>
void Density<T>::cutResolution(const T & min, const T & max){


	typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
	typedef typename Tarray::index index;
	index xdim = this->data3d.shape()[0];
	index ydim = this->data3d.shape()[1];
	index zdim = this->data3d.shape()[2];

	boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
	boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, zdim}};


        // where complex fft are stored
        Tarray_c  Afft(shapeft);
        // where real data is stored
        Tarray A(shape);

        for (int i = 0 ; i != xdim; ++i){
            for (int j = 0; j!=ydim; ++j){
                for (int k = 0; k != zdim; ++k){
                     A[i][j][k] = this->data3d[i][j][k] ;
                }
            }
        }

        // get dens3d ft
        metafft::FFT(A,Afft,-1);

        int n = 2;
        Afft[0][0][0] *=0.0;

        for (int i = -xdim/2 ; i != xdim/2; ++i){
            int ii = ((i<0)? i+xdim : i );
            for (int j = -ydim/2; j != ydim/2; ++j){
                int jj = ((j<0)? j+ydim : j);
                for (int k = -zdim/2;  k != zdim/2; ++k){
                    int kk = ((k<0)? k+zdim : k);
                    if (i != 0) {// just to make sure we don't divide by 0 when comp res.
                        Tvec v(3); v(1) = i; v(2) = j; v(3) = k;
                        T res = 1.0/std::sqrt(ublas::inner_prod(ublas::prod(v,this->imetrictensor),v));

                        // butterworth band filter in frequency (momentum) space
                        T w = (
                              (1.0/(1.0 + std::pow(min/res,2*n))) +
                              (1.0/(1.0 + std::pow(res/max,6*n)))
                              ) - 1.0;

                        Afft[ii][jj][kk] *= w;

                      }

                  }
              }
          }

        metafft::FFT(this->data3d,Afft,1);

        for (size_t i = 0 ; i!=this->data3d.num_elements(); ++i) {
            this->data3d.data()[i] /= this->data3d.num_elements();
          }

}




template <typename T>
void Density<T>::nlogtrans(){
    // check whether all the values are positive
    
    bool posi = true;
    
    size_t i = 0;
    
    while ( posi==true && (i < this->data3d.num_elements()) ) {
        
        if ( this->data3d.data()[i] < 0.0) {posi = false;}
        
        i++;
    }
    
    
    if (posi==true) {
        
        for (size_t i = 0; i != this->data3d.num_elements();++i){
            
            this->data3d.data()[i] = std::log(std::max(this->data3d.data()[i], 1e-12)) ;
        }
    } else {
        std::cout << "# Data is not positive." << std::endl;
    }
}



template <typename T>
void Density<T>::blobs_periodic(Density<T> & laplacian, const T & threshold) {

  T min =  *std::min_element(laplacian.data3d.data(),laplacian.data3d.data()+laplacian.data3d.num_elements()) ;


  std::cout <<"#\n# Laplacian analysis for a periodic system." << std::endl;
  std::cout << "# Volumetric data (rho): " << this->getName() << std::endl;
  std::cout << "# Laplacian L[rho]: " << laplacian.getName() << std::endl;
  std::cout << "# Considering zones with L[rho] < " << threshold << " * min(L[rho]) i.e. "<< min  << "." << std::endl;
  min *= threshold;
  std::cout << "# Considering zones with L[rho] < " << min << " [rho unit]/A^3" << std::endl;
  std::cout << "#\n#\n#" << std::endl;


  std::vector<std::vector<int>> lindeces;
  std::vector<double> vdensities, vlaplacian ;

  typedef typename Tarray::index index;
  index xdim = this->data3d.shape()[0];
  index ydim = this->data3d.shape()[1];
  index zdim = this->data3d.shape()[2];


  //std::cout << "Latice has: "<< xdim <<"x"<<ydim<<"x"<<zdim << std::endl;



  // (1) filter the vertices to be included based on the supplied laplacian threshold
  // and based on the value of density as some crystal densities contain negative values
  for (int i = 0 ; i != xdim; ++i){
      for (int j = 0; j!=ydim; ++j){
          for (int k = 0; k != zdim; ++k){
              Tvec v(3,0);
              v(0) = i;
              v(1) = j;
              v(2) = k;

              v = ublas::prod(this->rotm,v) + this->orig;

              double td  = this->data3d[i][j][k];
              double th = laplacian.get(v) ;

              if (th < min and td >= 0) {
                  std::vector<int> vv = {i, j, k};
                  lindeces.push_back(vv);
                  vdensities.push_back(td);
                  vlaplacian.push_back(th);
                }
            }
        }
    }

    
    // (2) build the adjacency matrix
    typedef std::pair<int,int> Edge;
    std::vector<Edge> used_by_v;
    
    
    for (size_t i = 0 ; i != lindeces.size();++i) {
        for (size_t j = i+1; j != lindeces.size();++j ) {

            int tmp0 = lindeces[i][0] - lindeces[j][0];
            int tmp1 = lindeces[i][1] - lindeces[j][1];
            int tmp2 = lindeces[i][2] - lindeces[j][2];
            
            // wrapping back
            if (std::abs(tmp0) == xdim-1) {tmp0=1;}
            if (std::abs(tmp1) == ydim-1) {tmp1=1;}
            if (std::abs(tmp2) == zdim-1) {tmp2=1;}

            // generating directed graph edges.
            if ( std::abs(tmp0) + std::abs(tmp1) + std::abs(tmp2) <= 3 ) {
                if (vlaplacian[i] <= vlaplacian[j])
                {
                    used_by_v.push_back(Edge(j,i));
                } else {
                    used_by_v.push_back(Edge(i,j));
                }
            }
            
            
        }
    }
    
    

    // (3) Init directed and undirected graphs to help
    // separate graph in connected components and map critical pointss
    typedef boost::adjacency_list
    <
        boost::vecS, boost::vecS, boost::undirectedS,
        boost::property< boost::vertex_color_t, int >
    > Graph;
    
    typedef boost::adjacency_list
    <
        boost::vecS, boost::vecS, boost::bidirectionalS,
        boost::property< boost::vertex_color_t, int >
    > DGraph;
    

    // init graphs
    Graph   g(used_by_v.data(), used_by_v.data() + used_by_v.size(), lindeces.size());
    DGraph dg(used_by_v.data(), used_by_v.data() + used_by_v.size(), lindeces.size());

    // checking wethers vertices are local minima
    // based on the number of outward and inwards edges
    std::vector<int> number_out_edges(0), number_in_edges(0);
    boost::graph_traits<DGraph>::vertex_iterator i, end;
    boost::graph_traits<DGraph>::in_edge_iterator ei, edge_end;
    boost::graph_traits<DGraph>::out_edge_iterator eei, eedge_end;
    
    for(boost::tie(i,end) = boost::vertices(dg); i != end; ++i) {
        boost::tie(eei,eedge_end) = boost::out_edges(*i, dg);
        number_out_edges.push_back(int(eedge_end - eei));
        
        boost::tie(ei,edge_end) = boost::in_edges(*i, dg);
        number_in_edges.push_back(int(edge_end - ei));
    }

    // separating graph in connected componnets
    std::vector<int> component(boost::num_vertices(g));
    int num_comp = boost::connected_components(
                                               g,
                                               boost::make_iterator_property_map(
                                                               component.begin(),
                                                               boost::get(boost::vertex_index, g)
                                                                                 )
                                               );
    
    

    std::string rootname = this->getName()+"-"+laplacian.getName();
    const std::string pdbfilecentr(rootname + "-blobs-centroid.pdb");
    const std::string pdbfilemean(rootname + "-blobs-meanposition.pdb");
    const std::string reportfile(rootname + "-blobs-data.dat");
    std::ofstream myfilep;
    std::ofstream myfilea;
    std::ofstream myfiled;
    
    
    //T factor = this->bulkc * 6.023 * 0.0001 * this->Vvox;
    T factor = this->bulkc * this->Vvox;
    myfilep.open(pdbfilecentr);
    myfilea.open(pdbfilemean);
//    myfiled.open(reportfile);
    
/*myfiled*/std::cout << "# Blob - index" << std::endl;
/*myfiled*/std::cout << "# Max[x,y,z] - coordinates of the point where density is maximum within the blob." << std::endl;
/*myfiled*/std::cout << "# Occ[#] - integrated occupancy of the blob. " << std::endl;
/*myfiled*/std::cout << "# V[A^3] - volume occupied by the blob." << std::endl;
/*myfiled*/std::cout << "# M[x,y,z] - average position (first moment) of the blob." << std::endl;
/*myfiled*/std::cout << "# Var[A^2] - average variance (second moment) of the blob." << std::endl;
/*myfiled*/std::cout << "# rave/max - ration between density at the average position and maximum value of density." << std::endl;
/*myfiled*/std::cout << boost::format("%-5s%-29s%9s%9s%-29s%9s%10s") %" Blob"%" Max[x,y,z]"%"Occ[#]"%"V[A^3]"%" M[x,y,z]"%"Var[A^2]"%" rave/max[%]" << std::endl;
    
    // index used for pdb output
    int atomindex = 1;
    
    // go over each graph component index
    for (int i = 0; i != num_comp ; ++i ){

        // integral over density
        T integral_density = 0;
        
        
        //T maxvalue = 0.0;
        Tvec maximum_density;
        
        // vector to store position (x,y,z)
        std::vector <Tvec> positions;
        
        // vector to store density
        std::vector <T> densities;
        std::vector<std::vector<int>> loclindeces (0);
        std::vector<int> loclmaxima (0);
        
        T maxvalue = -99999.0;
        size_t maxindex;
        
        for (size_t j = 0; j != component.size();++j) {
            if (component[j]==i) {
                if (number_out_edges[j]==0) {
                    loclmaxima.push_back(1);
                    
                    if (vdensities[j]>maxvalue){
                        maxvalue = vdensities[j];
                        maxindex = loclmaxima.size()-1;
                    }
                    
                }
                else{
                    loclmaxima.push_back(0);
                }
                loclindeces.push_back(lindeces[j]);
                densities.push_back(factor * vdensities[j]);
                integral_density += factor * vdensities[j];
            }
        }
        
        
        int num_local_maxima = std::accumulate(loclmaxima.begin(),loclmaxima.end(),int(0));
        
        //std::cout << "# local maxima: " << num_local_maxima << std::endl;
        
        if (num_local_maxima == 0) {
            T tt =  - 999999;
            //size_t index;
            for (size_t j = 0; j != densities.size(); ++j){
                if (densities[j] > tt ){
                    tt = densities[j];
                    maxindex = j;
                }
            }
            
            loclmaxima[maxindex] = 1;
        }


        // wrapping indexes
        for (size_t j = 0; j!=loclindeces.size(); ++j ){
           // std::cout << loclindeces[j] << " vs " ;
            if (j != maxindex){
                int tmp0 = loclindeces[maxindex][0] - loclindeces[j][0];
                int tmp1 = loclindeces[maxindex][1] - loclindeces[j][1];
                int tmp2 = loclindeces[maxindex][2] - loclindeces[j][2];

                int temp = tmp0 + xdim;
                if ( std::abs(temp) < std::abs(tmp0) ) {loclindeces[j][0] -= xdim; tmp0=temp;}
                temp = tmp0 - xdim;
                if ( std::abs(temp) < std::abs(tmp0) ) {loclindeces[j][0] += xdim; tmp0=temp;}

                temp = tmp1 + ydim;
                if ( std::abs(temp) < std::abs(tmp1) ) {loclindeces[j][1] -= ydim; tmp1=temp;}
                temp = tmp1 - ydim;
                if ( std::abs(temp) < std::abs(tmp1) ) {loclindeces[j][1] += ydim; tmp1=temp;}
            
                temp = tmp2 + zdim;
                if ( std::abs(temp) < std::abs(tmp2) ) {loclindeces[j][2] -= zdim; tmp2=temp;}
                temp = tmp2 - zdim;
                if ( std::abs(temp) < std::abs(tmp2) ) {loclindeces[j][2] += zdim; tmp2=temp;}
            }
            
            //std::cout << loclindeces[j] << std::endl;
            
        }

        // with indexes being wrapped, compute wrapped coordinates
        for (size_t j = 0; j!=loclindeces.size(); ++j ){
            Tvec v(3,0); v(0) = loclindeces[j][0] ; v(1) = loclindeces[j][1] ; v(2) = loclindeces[j][2];
            positions.push_back(ublas::prod(this->rotm,v) + this->orig);
        }

        // average position, first moment

        Tvec centroid_density(3,0.0);

        for (size_t t = 0;  t != positions.size(); ++t ){
            centroid_density += ( positions[t] * densities[t] )  ;
        }
        centroid_density /= integral_density;
        
        
        // variance

        Tvec variance_density(3,0.0);

        for (size_t t = 0;  t != positions.size(); ++t ){

            Tvec tmp = positions[t] - centroid_density;
            
            tmp(0) = tmp(0) * tmp(0);
            tmp(1) = tmp(1) * tmp(1);
            tmp(2) = tmp(2) * tmp(2);
            
            variance_density = variance_density + ( tmp * densities[t] )  ;
        }
        
        variance_density /= integral_density;
        T variance_density_ave = ublas::sum (variance_density)/3.0;
        
        for (size_t t =0; t != loclmaxima.size();++t){
            //std::cout << this->getper(positions[t]) << " " << positions[t] << " " << loclindeces[t] << std::endl;
            if (loclmaxima[t]==1)
            {
                // table print
               // myfiled <<
                //std::cout << "-----------------" << std::endl;
                //std::cout << this->getper(centroid_density) << " vs " << this->getper(positions[t]) << std::endl;
                //std::cout << "---------x-------" << std::endl;
                std::cout <<
                  boost::format ("%5i (%8.3f %8.3f %8.3f) %8.3f %8.3f (%8.3f %8.3f %8.3f) %8.3f %8.3f")
                  % i % positions[t][0] % positions[t][1] % positions[t][2]
                  % integral_density % (loclindeces.size() * this->Vvox)
                  % centroid_density[0] % centroid_density[1] % centroid_density[2]
                  % (0.515158 * std::pow(-1.0 * this->bulkc * laplacian.getper(positions[t]), -0.4))
                  % ((this->getper(centroid_density) / (maxvalue)))
                << std::endl;
          
              // pdb print
                myfilep  <<
                   boost::format("%-6s%5i%4s%2s%3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s")
                   % "ATOM" % atomindex % this->species % " " % this->residue % "w" % atomindex
                   % positions[t][0] % positions[t][1] % positions[t][2]
                   % integral_density
                   //% (40.6753 * std::pow(-1.0 * this->bulkc * laplacian.getper(positions[maxindex]), -0.4))
                   % std::fabs(laplacian.getper(positions[maxindex]))
                   % this->species
                << std::endl;
                
                atomindex++;
                
                // myfilep << "TER" << std::endl;
                
            }
        }
        // pdb average print
        myfilea  << boost::format("%-6s%5i%4s%2s%3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s")
            %"ATOM"%(i+1)%this->species%" "%this->residue%"w"%(i+1)
            //%centroid_density[0]%centroid_density[1]%centroid_density[2]
            % positions[maxindex][0] % positions[maxindex][1] % positions[maxindex][2]
            //%integral_density% (variance_density_ave * 26.3189)
            %integral_density
            //% (40.6753 * std::pow(-1.0 * this->bulkc * laplacian.get(positions[maxindex]), -0.4))
            % std::fabs(laplacian.getper(positions[maxindex]))
            %this->species
        << std::endl;
        // myfilea << "TER" << std::endl;


    }
    
    myfilep.close();
    myfilea.close();
//    myfiled.close();
    
    std::cout << "# Found " << num_comp << " blobs and " << atomindex << " tightly bound modes.\n#\n#" << std::endl;
    std::cout << "# Centroid coordinates in pdb format:" << std::endl;
    std::cout << "# \t"<< pdbfilecentr << std::endl;
    std::cout << "# Average coordinates of each blob in pdb format:" << std::endl;
    std::cout << "# \t"<< pdbfilemean << std::endl;
    std::cout << "# Tabulated detailed data for each blob:" << std::endl;
    std::cout << "# \t"<< reportfile << std::endl;

}


template <typename T>
void Density<T>::blobs(Density<T> & laplacian, const T & threshold, const T & bulkc) {
    
    std::cout <<"#\n# Laplacian analysis." << std::endl;
    std::cout <<"# ------------------------" << std::endl;
    
    size_t count = 0;
    
    
    std::vector<std::vector<int>> lindeces;
    std::vector<double> vdensities, vlaplacian ;
    
    typedef typename Tarray::index index;
    index xdim = this->data3d.shape()[0];
    index ydim = this->data3d.shape()[1];
    index zdim = this->data3d.shape()[2];
    
    // (1) filter the vertices to be included based on the supplied threshold
    for (int i = 0 ; i != xdim; ++i){
        for (int j = 0; j!=ydim; ++j){
            for (int k = 0; k != zdim; ++k){
                Tvec v(3,0);
                v(0) = i;
                v(1) = j;
                v(2) = k;
                
                v = ublas::prod(this->rotm,v) + this->orig;
                
                double td  = this->data3d[i][j][k];
                double th = laplacian.get(v) ;
                if (th < threshold) {
                    
                    count++;
                    std::vector<int> vv = {i, j, k};
                    lindeces.push_back(vv);
                    vdensities.push_back(td);
                    vlaplacian.push_back(th);
                }
            }
        }
    }
    
    
    // (2) build the adjacency matrix
    typedef std::pair<int,int> Edge;
    std::vector<Edge> used_by_v;
    
    for (size_t i = 0 ; i != lindeces.size();++i) {
        for (size_t j = i+1; j != lindeces.size();++j ) {
            
            int tmp0 = lindeces[i][0] - lindeces[j][0];
            int tmp1 = lindeces[i][1] - lindeces[j][1];
            int tmp2 = lindeces[i][2] - lindeces[j][2];
            
            if ( (tmp0*tmp0 + tmp1*tmp1 + tmp2*tmp2) == 1 ) {
                
                double t  = this->data3d[lindeces[i][0] ] [ lindeces[i][1]] [ lindeces[i][2]];
                double u  = this->data3d[lindeces[j][0] ] [ lindeces[j][1]] [ lindeces[j][2]];
                
                if ( t >= u ) {
                    used_by_v.push_back(Edge(j,i));
                }
                else
                {
                    used_by_v.push_back(Edge(i,j));
                }
            }
        }
    }
    
    
    
    // std::cout << "# Vertices under threshold: " << count << " out of " << xdim * ydim * zdim << "."<< std::endl;
    
    // type of graph
    
    typedef boost::adjacency_list
    <
    boost::vecS, boost::vecS, boost::undirectedS,
    boost::property< boost::vertex_color_t, int >
    > Graph;
    
    typedef boost::adjacency_list
    <
    boost::vecS, boost::vecS, boost::bidirectionalS,
    boost::property< boost::vertex_color_t, int >
    > DGraph;
    
    //    typedef typename boost::graph_traits<DGraph>::vertex_descriptor DVertex;
    // init graphs
    Graph g(used_by_v.data(), used_by_v.data() + used_by_v.size(), lindeces.size());
    DGraph dg(used_by_v.data(), used_by_v.data() + used_by_v.size(), lindeces.size());
    
    
    
    const std::string filenameout("blobs.graphml");
    boost::dynamic_properties dp;
    
    boost::graph_traits<Graph>::vertex_iterator v, v_end;
    for (boost::tie(v,v_end) = vertices(dg); v != v_end; ++v)
        boost::put(boost::vertex_color_t(), dg, *v, *v);
    /*
     graph_traits<Graph>::edge_iterator e, e_end;
     for (tie(e,e_end) = edges(g); e != e_end; ++e)
     put(edge_weight_t(), g, *e, 3);
     */
    
    dp.property("name", boost::get(boost::vertex_color_t(), dg));
    /*
     dp.property("weight", get(edge_weight_t(), g));
     */
    
    //        myfile.open(filenameout);
    //        boost::write_graphml(myfile, dg, dp, true);
    //        myfile.close ();
    
    std::vector<int> component(boost::num_vertices(g));
    
    int num_comp = boost::connected_components(
                                               g,
                                               boost::make_iterator_property_map(
                                                                                 component.begin(),
                                                                                 boost::get(boost::vertex_index, g))
                                               );
    
    
    // std::cout << "# number of graph components (subgraphs): " << num_comp << std::endl;
    
    
    std::string rootname = this->getName()+"-"+laplacian.getName();
    
    //rootname += laplacian.getName();
    
    //std::cout << "# >> rootname: " << rootname << std::endl;
    
    
    
    //const std::string pdbfile(rootname + "-blobs.pdb");
    const std::string pdbfilecentr(rootname + "-blobs-centroid.pdb");
    const std::string pdbfilemean(rootname + "-blobs-meanposition.pdb");
    const std::string reportfile(rootname + "-blobs-data.dat");
    
    
    //std::cout << "# "<< pdbfile << "\t contains some information." << std::endl;
    std::cout << "# << "<< pdbfilecentr << ": centroid coordinates in pdb format." << std::endl;
    std::cout << "# << "<< pdbfilemean << ": average coordinates of each drop in pdb format." << std::endl;
    std::cout << "# << "<< reportfile << ": tabulated detailed data for each drop." << std::endl;
    
    
   // std::ofstream myfile;
    std::ofstream myfilep;
    std::ofstream myfilea;
    std::ofstream myfiled;
    
    
    
    //double voxel_volume = this->rotm(0,0) * this->rotm(1,1) * this->rotm(2,2);
    //double factor = this->bulkc * 6.023 * 0.0001 * this->Vvox;
    double factor = this->bulkc * this->Vvox;
    //myfile.open(pdbfile);
    myfilep.open(pdbfilecentr);
    myfilea.open(pdbfilemean);
    myfiled.open(reportfile);
    
    myfiled << "# Blob - index" << std::endl;
    myfiled << "# Max[x,y,z] - coordinates of the point where density is maximum within the blob." << std::endl;
    myfiled << "# Occ[#] - integrated occupancy of the blob. " << std::endl;
    myfiled << "# V[A^3] - volume occupied by the blob." << std::endl;
    myfiled << "# M[x,y,z] - average position (first moment) of the blob." << std::endl;
    myfiled << "# Var[A^2] - average variance (second moment) of the blob." << std::endl;
    myfiled << "# rave/max - ration between density at the average position and maximum value of density." << std::endl;
    myfiled << boost::format("%-5s%-29s%9s%9s%-29s%9s%10s") %" Blob"%" Max[x,y,z]"%"Occ[#]"%"V[A^3]"%" M[x,y,z]"%"Var[A^2]"%" rave/max[%]" << std::endl;
    
    // go over each graph component index
    for (int i = 0; i != num_comp ; ++i ){
        
        // integral over density
        double integral_density = 0;
        
        // centroid coordinates
        Tvec centroid_density(3,0.0);
        Tvec variance_density(3,0.0);
        
        T maxvalue = 0.0;
        Tvec maximum_density;
        
        
        //       T sumweight = 0.0;
        
        int nvoxels = 0;
        
        
        // vector to store position (x,y,z)
        
        std::vector <Tvec> positions;
        
        // vector to store density
        
        std::vector <double> densities;
        
        // go over each vertex and check whether is in the component
        for (int j = 0; j != component.size();++j) {
            
            
            if (component[j]==i) {
                
                Tvec v(3,0); v(0) = lindeces[j][0] ; v(1) = lindeces[j][1] ; v(2) = lindeces[j][2];
                
                v = ublas::prod(this->rotm,v) + this->orig;
                
                T weight = factor*vdensities[j] ;
                
                
                positions.push_back(v);
                densities.push_back(weight);
                
                
                /*
                myfile <<
                boost::format("ATOM  %5i  XXX   Y Z %3i    %8.3f%8.3f%8.3f  1.00%6.2f")
                %j %i %v(0) %v(1) %v(2) %(weight*100)
                << std::endl;
                */
                
                integral_density += weight;
                
                // keeping record of the max value so far.
                if (maxvalue < vdensities[j]) {
                    maxvalue=vdensities[j] ;
                    maximum_density = v;
                }
                
                // counting the number of voxels
                ++nvoxels;
            }
            
            
            
        }
        
        
        
        for (size_t t = 0;  t != positions.size(); ++t ){
            
            centroid_density = centroid_density + ( positions[t] * densities[t] )  ;
            
        }
        
        centroid_density /= integral_density;
        
        
        
        for (size_t t = 0;  t != positions.size(); ++t ){
            
            
            Tvec tmp = positions[t] - centroid_density;
            
            tmp(0) = tmp(0) * tmp(0);
            tmp(1) = tmp(1) * tmp(1);
            tmp(2) = tmp(2) * tmp(2);
            
            variance_density = variance_density + ( tmp * densities[t] )  ;
            
        }
        
        variance_density /= integral_density;
        variance_density(0) = std::sqrt( variance_density(0) );
        variance_density(1) = std::sqrt( variance_density(1) );
        variance_density(2) = std::sqrt( variance_density(2) );
        double variance_density_ave = ublas::norm_2(variance_density)* std::sqrt(1.0/3.0);
        
        
        std::cout << std::fixed << std::setprecision(3) ;
        //myfile << "TER" << std::endl;
        
        if (integral_density > 0.001) {//outputing only useful information,
            //            std::cout << "# Blob_" << i << " ";
            //            std::cout << centroid_density << " +/- " <<  variance_density_ave  << " ";
            //            std::cout << " ->(%_of_max) " << (this->get(centroid_density) / maxvalue) *100.0  << " ";
            //            std::cout << " max_@ " << maximum_density ;
            //            std::cout << " Occupation[#] = " << integral_density  << " ";
            //            std::cout << " Blob_Volume[A^3] = " << nvoxels * voxel_volume  <<  std::endl;
            
            
            myfiled << boost::format ("%5i (%8.3f %8.3f %8.3f) %8.3f %8.3f (%8.3f %8.3f %8.3f) %8.3f %8.3f") %i %centroid_density[0] %centroid_density[1] %centroid_density[2] %integral_density %(nvoxels * this->Vvox) %maximum_density[0] %maximum_density[1] %maximum_density[2] %variance_density_ave %( (this->get(centroid_density) / maxvalue) *100.0) << std::endl;
            
            myfilep <<
            boost::format("ATOM  %5i  XXX   Y Z %3i    %8.3f%8.3f%8.3f  1.00%6.2f")
            %i %i %maximum_density(0) %maximum_density(1) %maximum_density(2) %(integral_density)
            << std::endl;
            myfilep << "TER" << std::endl;
            
            
            myfilea <<
            boost::format("ATOM  %5i  XXX   Y Z %3i    %8.3f%8.3f%8.3f  1.00%6.2f")
            %i %i %centroid_density(0) %centroid_density(1) %centroid_density(2) %(integral_density)
            << std::endl;
            myfilea << "TER" << std::endl;
        }
        
        
    }
    
    
    //myfile.close();
    myfilep.close();
    myfilea.close();
    myfiled.close();
    
    std::cout << "#\n# -------------------------\n#" << std::endl;
    
    
    // Write out the outgoing edges
    //typename boost::graph_traits<DGraph>::out_edge_iterator out_i, out_end;
    //typename boost::graph_traits<DGraph>::edge_descriptor e;
    
    //            boost::graph_traits<DGraph>::vertex_iterator v, v_end;
    /*            for (boost::tie(v,v_end) = boost::vertices(dg); v != v_end; ++v)
     {
     std::cout << "(v"<< *v << ")\t" << boost::out_degree(*v,dg) << "\t " << boost::in_degree(*v,dg) << std::endl;
     }
     */
    /*
     for (boost::tie(v,v_end) = boost::vertices(dg); v != v_end; ++v)
     {
     if (boost::out_degree(*v,dg)==0){
     std::cout << "(local max v "<< *v << " )\t" << boost::out_degree(*v,dg) << "\t " << boost::in_degree(*v,dg) << std::endl;
     }
     }
     */
    
}

template <typename T>
void Density<T>::writedxfile(const std::string &filenameout){

  std::vector<std::string> stokens;
  stokens =  boost::copy_range<std::vector<std::string>>(filenameout |boost::adaptors::tokenized( std::string("\\."), -1 ) );
  std::string ext;
  ext = stokens[stokens.size()-1];


  if (ext=="dx") {

      std::cout << "# >> Writing density to " << filenameout  << "." << std::endl;

      std::ofstream myfile;
      myfile.open(filenameout);
      //myfile.precision(5);
      //myfile << std::scientific;

      myfile << boost::format("# generated by moft.metaTwist\n") ;
      myfile << boost::format("object 1 class gridpositions counts %8.3f %8.3f %8.3f\n") % this->data3d.shape()[0] %  this->data3d.shape()[1]  % this->data3d.shape()[2]  ;
      myfile << boost::format("origin %8.3f %8.3f %8.3f\n") %  this->orig(0)   %  this->orig(1)  % this->orig(2);
      myfile << boost::format("delta  %8.3f %8.3f %8.3f\n")  %  this->rotm(0,0) % this->rotm(1,0) % this->rotm(2,0);
      myfile << boost::format("delta  %8.3f %8.3f %8.3f\n")  %  this->rotm(0,1) % this->rotm(1,1) % this->rotm(2,1);
      myfile << boost::format("delta  %8.3f %8.3f %8.3f\n")  %  this->rotm(0,2) % this->rotm(1,2) % this->rotm(2,2);
      myfile << boost::format("object 2 class gridconnections counts %8i %8i %8i \n")% this->data3d.shape()[0] %  this->data3d.shape()[1]  % this->data3d.shape()[2] ;
      myfile << boost::format("object 3 class array type double rank 0 items %8i data follows\n") % (this->data3d.shape()[0] *  this->data3d.shape()[1] *  this->data3d.shape()[2]);

      for (size_t i = 0; i != this->data3d.num_elements();++i){

          myfile << boost::format("%12.4E ") % this->data3d.data()[i];

          if ( (i !=0) && ((i+1)%3==0)) myfile << "\n";

        }

      myfile << "\nobject \"Untitled\" class field" ;
    }

  else if (ext=="ccp4" || ext=="mrc") {

      this->writeccp4(filenameout);

    } else {

      std::cout << "# Extension " << ext << " not recognized." << std::endl;
      std::string newfilename;
      newfilename =  filenameout+std::string(".ccp4");
      std::cout << "# Writing to " <<  newfilename << "." << std::endl;
      this->writeccp4(newfilename);

    }

}



template <typename T>
void Density<T>::writeccp4(const std::string &filenameout){
      std::cout << "# Writing volumetric data to " << filenameout << "." << std::endl;
      util::files::writeccp4(filenameout, this->orig, this->a, this->b, this->c, this->data3d);
}

template <typename T>

void Density<T>::getQDensity(Density<T> & QDensity){

    for (size_t i = 0; i != this->data3d.num_elements(); ++i) {
        QDensity.data3d.data()[i] = this->data3d.data()[i] * this->partialcharge * this->bulkc;
    }

}

template <typename T>
void Density<T>::getLaplacian(Density<T> & LDensity){
    
    
    // Check if grid sizes and coordinate system match.
    
    
    std::cout << "# >> Computing numerical Laplacian of input density." << std::endl;
    
    if (LDensity.getorig() == this->getorig()
        && LDensity.getrotm() == this->getrotm()
        && LDensity.getGridSize() == this->getGridSize()
        )
    {
        //std::cout << "1" << std::endl;
        
        int isize = this->data3d.shape()[0], jsize = this->data3d.shape()[1], ksize = this->data3d.shape()[2];
        
        
        T factor = this->bulkc * 6.0/std::pow(this->rotm(1,1),2);
        
        //std::cout << "factor = " << factor << std::endl;
        
        for (int i = 0; i != isize; ++i){
            int im = i-1, ip=i+1;
            
            for (int j = 0; j != jsize; ++j) {
                int jm = j-1, jp=j+1;
                
                for (int k = 0; k != ksize; ++k) {
                    int km = k -1, kp = k+1;
                    
                    
                    size_t count = 0;
                    
                    
                    T depsilon = 0;
                    T d0 = this->data3d[i][j][k];
                    
                    if (im>=0)    {count++; depsilon += this->data3d[im][j ][k ]; }
                    if (jm>=0)    {count++; depsilon += this->data3d[i ][jm][k ]; }
                    if (km>=0)    {count++; depsilon += this->data3d[i ][j ][km]; }
                    if (ip<isize) {count++; depsilon += this->data3d[ip][j ][k ];}
                    if (jp<jsize) {count++; depsilon += this->data3d[i ][jp][k ];}
                    if (kp<ksize) {count++; depsilon += this->data3d[i ][j ][kp];}
                    
                    
                    depsilon /= T(count);
                    
                    //
                    //                    if (std::abs(depsilon-d0) > T(1000.0) ){
                    //
                    //                        std::cout << i << " " << j << " " << k << " " << depsilon << " - " << d0 << " = " <<(depsilon - d0)*factor << std::endl;
                    //                    }
                    
                    
                    LDensity.set(i,j,k,factor*(depsilon  - d0));
                    
                    
                }
            }
            
            
            
        }
        
        
    } else {
        
        std::cout << "# Laplacian density does not have the same size as the density object. Skipping.`" << std::endl;
        
    }
    
    
    // go over the interior grid points
    
    // go over the points at the grid margins
    
}


template <typename T>
void Density<T>::getNLaplacian(Density<T> & LDensity){
    
    // Check if grid sizes and coordinate system match.
    
    std::cout << "# >> Computing numerical Normalized Laplacian of input density." << std::endl;
    
    if (LDensity.getorig() == this->getorig()
        && LDensity.getrotm() == this->getrotm()
        && LDensity.getGridSize() == this->getGridSize()
        )
    {
        //std::cout << "1" << std::endl;
        
        int isize = this->data3d.shape()[0], jsize = this->data3d.shape()[1], ksize = this->data3d.shape()[2];
        
        
        //T factor = 6.0/std::pow(this->rotm(1,1),2);
        
        //std::cout << "factor = " << factor << std::endl;
        
        for (int i = 0; i != isize; ++i){
            int im = i-1, ip=i+1;
            
            for (int j = 0; j != jsize; ++j) {
                int jm = j-1, jp=j+1;
                
                for (int k = 0; k != ksize; ++k) {
                    int km = k -1, kp = k+1;
                    
                    
                    size_t count = 0;
                    
                    
                    T depsilon = 0.0;
                    T d0 = this->data3d[i][j][k];
                    
                    if (im>=0)    {count++; depsilon += this->data3d[im][j ][k ];}
                    if (jm>=0)    {count++; depsilon += this->data3d[i ][jm][k ];}
                    if (km>=0)    {count++; depsilon += this->data3d[i ][j ][km];}
                    if (ip<isize) {count++; depsilon += this->data3d[ip][j ][k ];}
                    if (jp<jsize) {count++; depsilon += this->data3d[i ][jp][k ];}
                    if (kp<ksize) {count++; depsilon += this->data3d[i ][j ][kp];}
                    
                    // Average density in neighbourhood
                    depsilon /= T(count);
                    // Numerical safety.
                    depsilon = std::max(depsilon,1e-3);
                    
                    T tmp = /*factor * */( 1.0 - ( d0/depsilon ));
                    
                    //if ( std::abs(tmp)  > 50) {
                    //    std::cout << depsilon << " vs " << d0 << " = " << depsilon/std::max(d0,1e-6) << std::endl;
                    //}
                    
                    //T tmp = depsilon/std::max(d0,1e-3);
                    
                    LDensity.set(i,j,k,tmp);
                    
                }
            }
            
            
            
        }
        
        
    } else {
        
        std::cout << "Laplacian density does not have the same size as the density object. Skipping.`" << std::endl;
        
    }
    
    
    

    
    
    // go over the interior grid points
    
    // go over the points at the grid margins
    
}
template <typename T>
void Density<T>::getNpLaplacian(Density<T> & LDensity) {
    
    // Check if grid sizes and coordinate system match.
    
    
    std::cout << "# >> Computing numerical Normalized Laplacian of input density." << std::endl;
    
    if (LDensity.getorig() == this->getorig()
        && LDensity.getrotm() == this->getrotm()
        && LDensity.getGridSize() == this->getGridSize()
        ) {
        //std::cout << "1" << std::endl;
        
        int isize = this->data3d.shape()[0], jsize = this->data3d.shape()[1], ksize = this->data3d.shape()[2];
        
        
 //       T factor = 6.0 / std::pow(this->rotm(1, 1), 2);
        
        //std::cout << "factor = " << factor << std::endl;
        
        for (int i = 0; i != isize; ++i) {
            int im = i - 1, ip = i + 1;
            
            for (int j = 0; j != jsize; ++j) {
                int jm = j - 1, jp = j + 1;
                
                for (int k = 0; k != ksize; ++k) {
                    int km = k - 1, kp = k + 1;
                    
                    
                    size_t count = 0;
                    
                    
                    T depsilon = 0.0;
                    T d0 = this->data3d[i][j][k];
                    
                    if (im >= 0) {
                        count++;
                        depsilon += this->data3d[im][j][k];
                    }
                    if (jm >= 0) {
                        count++;
                        depsilon += this->data3d[i][jm][k];
                    }
                    if (km >= 0) {
                        count++;
                        depsilon += this->data3d[i][j][km];
                    }
                    if (ip < isize) {
                        count++;
                        depsilon += this->data3d[ip][j][k];
                    }
                    if (jp < jsize) {
                        count++;
                        depsilon += this->data3d[i][jp][k];
                    }
                    if (kp < ksize) {
                        count++;
                        depsilon += this->data3d[i][j][kp];
                    }
                    
                    // Average density in neighbourhood
                    depsilon /= T(count);
                    // Numerical safety.
                    depsilon = std::max(depsilon, 1e-3);
                    
                    T tmp =  std::log( std::max(d0/depsilon,1e-8) );
                    
                    //if ( std::abs(tmp)  > 50) {
                    //    std::cout << depsilon << " vs " << d0 << " = " << depsilon/std::max(d0,1e-6) << std::endl;
                    //}
                    
                    //T tmp = depsilon/std::max(d0,1e-3);
                    
                    LDensity.set(i, j, k, tmp);
                    
                    
                }
            }
            
            
        }
        
        
    } else {
        
        std::cout << "Laplacian density does not have the same size as the density object. Skipping.`" << std::endl;
        
    }
}






#endif /* DENSITY_HPP_ */
