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


/*

  This is metaTWIST: a simple low level driver for most of the functionalities available
  in MoFT.

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <metatwist/boosthead.hpp>
#include <utilities/utilities.hpp>
#include <utilities/rotations.hpp>
#include <utilities/interpolation.hpp>
#include <metatwist/Density.hpp>


typedef boost::random::mt19937 Tmt;

int main(int ac, char* av[]) {
    
    typedef ublas::matrix<double> TMat;
    typedef ublas::vector<double> TVec;
    
    std::string filename,maptype,filenameout,/*filenamehess,*/ stringfile;
    double rmax,zmax,xmax,ymax;
    double bin;
    double utrate; // untwisting rate
    double bulkdens;
    double charge;
    double sigma;
    double threshold;
    size_t convtype;
    std::vector<double> com(3,0);
    std::vector< std::string > dxfilenames, filenamehess;

    std::vector<std::string> chemicalspecies;

    std::vector<double> resminmax = {1.0,10.0} ;
    
    std::string date(__DATE__);
    date +=" @ ";
    date +=__TIME__;

    unsigned int width = 100;
    po::options_description options
        ("\n\nThis is metaTWIST:  a simple low level driver for most of the functionality available in MoFT. \n\nCompiled on: "+date+".\n\nOptions",
         width,
         width/4);
    
    options.add_options()
    ("help,h",  "Produces help message.")
    //("dx",      po::value<std::string>(&filename)->required(), 	"DX file")
    ("dx",      po::value< std::vector<std::string> >()->required()->multitoken(),
                "Input density file(s): *.dx(gz,bz2)|*.ccp4."
    )
    ("map",     po::value<std::string>(&maptype)/*->required()*/,
     
     "Mapping type:\
     \n~ cylindrical (1D): cylindrical RDF along z-axis.\
     \n~ twist (2D): twisted helical map along z-axis.\
     \n~ untwist (2D): untwisted helical map along z-axis.\
     \n~ spherical (1D): spherical RDF.\
     \n~ projxyz: (1D) project 3D-map on x,y,z axes.\
     \n~ excess: excess number of particles.\
     \n~ blobs: Laplacian blob analysis.\
     \n~ blobsper: Laplacian blob analysis on a periodic 3D-map.\
     \n~ rhoel: Electron density using atomic form factors.\
     \n~ rhoelreal: Electron density using atomic densities.\
     \n~ charge: Partial charge.\
     \n~ cutresol: Cut 3D-map resolution range.\
     "
     )
    ("bin",     po::value<double>(&bin)->default_value(0.1, "0.1"),
     "Bin size for re-sampling."
     )
    ("rmax",    po::value<double>(&rmax)->default_value(60.0),
     "Extent along the rho direction (cylindrical, spherical RDFs)."
     )
    ("zmax",    po::value<double>(&zmax)->default_value(70.0),
     "Extent in the z direction."
     )
    ("ymax",    po::value<double>(&ymax)->default_value(20.0),
     "Extent in the y direction."
     )
    ("xmax",    po::value<double>(&xmax)->default_value(20.0),
     "Extent in the x direction."
     )
    ("utrate",  po::value<double>(&utrate)->default_value(0.18587, "0.18587"),
     "Untwisting rate: \
     \n0.18587 rad/Ang - BDNA(default)\
     \n0.16870 rad/Ang - TDNA\
     \n0.25590 rad/Ang - ARNA."
     )
    ("com",     po::value< std::vector<double> > (&com)->multitoken(),
     "COM coordinates."
     )
    ("resolution", po::value< std::vector<double> > (&resminmax)->multitoken(),
     "Min and max resolution thresholds (default: 1.0 and 10.0 Ang)."
      )
    ("bulkdens",po::value<double>(&bulkdens)->default_value(1.0),
     "Bulk density (A^-3)."
     )
    ("species",po::value<std::vector<std::string>>(&chemicalspecies)->multitoken()->required(),
     "Chemical species: atom, e.g. N, or atom & residue, e.g. O WAT )."
     )
     ("charge",po::value<double>(&charge),
       "Partial charge.)."
     )
    ("sigma",   po::value<double>(&sigma)->default_value(0.0),
     "Convolution sigma."
     )
    ("threshold", po::value<double>(&threshold)->default_value(0.0),
     "Laplacian threshold."
     )
    ("convolve",po::value<size_t>(&convtype)->default_value(4),
     "Convolution type\n - (1) Gaussian\n - (2) box\n - (3) sinc\n - (4) Laplacian of Gaussian."
     )
    ("odx",     po::value<std::string>(&filenameout),
     "Output density file (*.dx|*.ccp4)."
     )
    //("ldx",     po::value<std::string>(&filenamehess),
    ("ldx",     po::value< std::vector<std::string> >()->multitoken(),
     "Input Laplacian file (*.dx|*.ccp4)."
     )
    ("nlog",
     "Take the negative natural logarithm of the input density."
     )
    ("laplacian",
     "Compute Laplacian, L[rho], of the input density using finite difference."
     )
//    ("nlaplacian",
//     "Compute the normalized Laplacian, L[rho(x,y,z)]/rho_epsilon(x,y,z) of the input density using finite difference."
//     )
//    ("nplaplacian",
//     "Compute the normalized Laplacian, log2[rho(x,y,z)/rho_epsilon(x,y,z)] of the input density using finite difference."
//     )
//    ("worm",  po::value<std::string>(&stringfile),
//     " File containing the string coordinates."
//     )
    ("average","Average volumetric data, in case multiple datasets have been loaded. Otherwise, data will be accumulated.")
    ;
    // parsing input:
    po::variables_map vm;
    po::parsed_options parsed_options = po::parse_command_line(ac, av, options, po::command_line_style::unix_style ^ po::command_line_style::allow_short);
    //po::store(po::parse_command_line(ac, av, options, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::store(parsed_options, vm);
    
    
    if (vm.count("help")) {
        std::cout << options << "\n";
        return 1;
    }
    
    if (ac > 1) {
        // add date, user, machine name, pwd
        std::cout << "#\n#\n# metaTwist (compiled on " <<date << ")\n#\n#" <<std::endl;
        std::cout << "# Reading from input: \n# metatwist " ;
        for (int i = 1; i != ac; ++i) std::cout << av[i] << " ";
        std::cout << std::endl;
        std::cout << "#\n#" << std::endl;
    }
    
    // make sure the parsed file names are accumulated
    for (const po::option& o : parsed_options.options) {
        if (o.string_key == "dx") {
            dxfilenames = o.value;
        }
        
        if (o.string_key == "ldx"){
            filenamehess = o.value;
        }
    }
    
    // input checks:
    try{
        po::notify(vm);
        if (rmax <= 0.0 || zmax <=0 ){
            std::cout << options << "\n";
        }
    }
    catch(std::exception& e)
    {
        std::cout << options << "\n";
        std::cout << e.what() << "\n";
        return 1;
    }
    
    std::vector< TVec > allvcom;
    TVec vcom (3,0);
    // std::cout << com << std::endl;
    
    
    // parsing origin/center of mass information
    if (com.size()%3==0){
        //std::cout << "# Reading center of mass info." << std::endl;
        size_t t = com.size() / 3;
        size_t cnt = 0;
        allvcom.resize(t);
        for (size_t i = 0; i < t; ++i){
            allvcom[i].resize(3);
            allvcom[i](0) = com[cnt]; cnt++;
            allvcom[i](1) = com[cnt]; cnt++;
            allvcom[i](2) = com[cnt]; cnt++;
            //std::cout << "# "<< allvcom[i] << std::endl;
        }
        //std::cout << "# Done reading center of mass info." << std::endl;
    }
    
    //	if (com.size() == 3){
    //		for(size_t i = 0 ; i != vcom.size(); ++i)vcom(i)=com[i];
    //	}
    
    vcom = allvcom[0];
    // std::cout << "# com " << vcom << std::endl;
    
    typedef Density<double> Tdens;
    
    
    // Tdens   dens3d(filename);
    Tdens   dens3d(dxfilenames,chemicalspecies);
    
    
    
    // condition moving fwd to having read the density file
    
    if (!dens3d.getreadany()){
        
        std::cout << "# [x] Input density file could not been read. Exiting..." << std::endl;
        
        return 1;
        
    }
    
    dens3d.setBulkC(bulkdens);
    dens3d.printDetails();
    
    
    //dens3d.setspecies(chemicalspecies);
    
    if (vm.count("average")) {
        std::cout << "# Averaging loaded data: (d[1]+d[2]+...+d[n])/n." << std::endl;
        dens3d.average();

     }


    
    if (sigma > 0.0) {
        std::cout << "# Convolving density map." << std::endl;
        dens3d.convolutex(sigma,convtype);
        dens3d.writedxfile(filenameout);
    }
    
    if (vm.count("laplacian")) {
        Tdens lapldens3d(dens3d);
        dens3d.getLaplacian(lapldens3d);
        lapldens3d.writedxfile(std::string("laplacian-")+filenameout);
    }
    
    if (vm.count("nlaplacian")) {
        Tdens lapldens3d(dens3d);
        dens3d.getNLaplacian(lapldens3d);
        lapldens3d.writedxfile(std::string("nlaplacian-")+filenameout);
    }
    
    if (vm.count("nplaplacian")) {
        Tdens lapldens3d(dens3d);
        dens3d.getNpLaplacian(lapldens3d);
        lapldens3d.writedxfile(std::string("nplaplacian-")+filenameout);
    }
    
    if (vm.count("nlog")) {
        std::cout << "# Taking the negative logarithm." << std::endl;
        dens3d.nlogtrans();
        dens3d.writedxfile(std::string("nlog-")+filenameout);
    }


    if (vm.count("charge")) {
        std::cout << "# Computing charge density." << std::endl;
        dens3d.setcharge(charge);
        Tdens qdensity(dens3d);
        dens3d.getQDensity(qdensity);
        qdensity.writedxfile(std::string("qdens-")+filenameout);
    }

    if (vm.count("odx")
        and
        !( sigma>0.0 // aka convolve
          or
          vm.count("laplacian") or
          vm.count("nlaplacian") or
          vm.count("nplaplacian") or
          vm.count("nlog") or
          vm.count("map")
          )
        )
    {
        std::cout << "# Writing density to file." << "\n";
        dens3d.writedxfile(filenameout);
    }
    
    
    
    
    if (maptype=="cylindrical"){
        std::cout << "# 1D-cylindrical mapping:\n Remap the input 3D map in cylindrical coordinates where \nthe cylindrical axis is along the z-axis." << std::endl;
        
        
        size_t nz = 10, nr = 10, nt = 10;
        double z0 = -1.0 * zmax;
        double r0 = 0.0;
        double t0 = 0.0,tmax = 2.0*M_PI;
        
        double zinc =0.1 , tinc = M_PI/90.0;
        
        nt = size_t((tmax - t0)/tinc);
        nz = size_t((zmax - z0)/zinc);
        nr = size_t((rmax - r0)/bin);
        
        double r,z,t;
        
        TVec grid(nr,0.0),gridc(nr,0.0);
        
        for (size_t ir = 0; ir != nr; ++ir){
            r = r0 + ir * bin;
            for (size_t iz = 0; iz != nz; ++iz){
                z =  z0 + iz * zinc;
                for (size_t it = 0; it != nt; ++it){
                    t = t0 + it * tinc;
                    TVec v(3, 0.0);
                    v(0) = r*std::cos(t);
                    v(1) = r*std::sin(t);
                    v(2) = z;
                    grid[ir]  += dens3d.get(v+vcom);
                    gridc[ir] += 1.0;
                    
                }
            }
        }
        
        
        for (size_t i = 0; i != grid.size(); ++i){
            if (gridc[i]>0){
                std::cout << i*bin << " " << grid[i] / gridc[i] << std::endl;
            }
            
        }
        
    } else if (maptype=="twist")
    {
        
        std::cout << "# 2D-cylindrical twisted mapping." << std::endl;
        // number of bins
        
        size_t nx = size_t(2*xmax/bin);
        size_t ny = size_t(2*ymax/bin);
        size_t nz = size_t(2*zmax/bin);
        
        double x0 = -1.0*xmax;
        double y0 = -1.0*ymax;
        double z0 = -1.0*zmax;
        
        TMat grid(nx,ny,0),gridc(nx,ny,0);
        
        for (size_t ix = 0; ix != nx; ++ix){
            for (size_t iy = 0; iy != ny; ++iy){
                for (size_t iz = 0; iz != nz; ++iz){
                    
                    TVec v(3, 0.0);
                    v(0) = x0+ix*bin;
                    v(1) = y0+iy*bin;
                    v(2) = z0+iz*bin;
                    v+= vcom;
                    grid(ix,iy)  += dens3d.get(v);
                    gridc(ix,iy) += 1.0;
                    
                }
            }
        }
        
        
        
        for (size_t i = 0; i != grid.size1();++i){
            for (size_t j = 0; j != grid.size2();++j){
                std::cout << x0+i*bin << " " << y0+j*bin << " ";
                if(gridc(i,j)!=0.0){
                    std::cout << grid(i,j)/gridc(i,j) << std::endl;
                } else {
                    std::cout << 0.0 << std::endl;
                }
                
            }
            
            std::cout << std::endl;
        }
        
        
        
        
    } else if (maptype=="untwist"){
        
        
        std::cout << "# 2D-cylindrical untwisted mapping." << std::endl;
        // number of bins
        
        size_t nx = size_t(2*xmax/bin);
        size_t ny = size_t(2*ymax/bin);
        size_t nz = size_t(2*zmax/bin);
        
        double x0 = -1.0*xmax;
        double y0 = -1.0*ymax;
        double z0 = -1.0*zmax;
        
        TMat grid(nx,ny,0),gridc(nx,ny,0);
        
        
        TVec zax(3,0);
        zax(2) = 1.0;
        TMat rot;
        
        for (size_t iz = 0; iz != nz; ++iz){
            
            double theta = utrate* (z0+iz*bin) ;
            rotations::rotmatrixaxis(theta,zax,rot);
            
            for (size_t ix = 0; ix != nx; ++ix){
                for (size_t iy = 0; iy != ny; ++iy){
                    
                    
                    TVec v(3, 0.0);
                    v(0) = x0+ix*bin;
                    v(1) = y0+iy*bin;
                    v(2) = z0+iz*bin;
                    grid(ix,iy)  += dens3d.get(ublas::prod(rot, v + vcom));
                    gridc(ix,iy) += 1.0;
                    
                }
            }
        }
        
        
        
        for (size_t i = 0; i != grid.size1();++i){
            for (size_t j = 0; j != grid.size2();++j){
                std::cout << x0+i*bin << " " << y0+j*bin << " ";
                if(gridc(i,j)!=0.0){
                    std::cout << grid(i,j)/gridc(i,j) << std::endl;
                } else {
                    std::cout << 0.0 << std::endl;
                }
                
            }
            
            std::cout << std::endl;
        }
        
    }
    
    else if (maptype=="spherical"){
        
        
        std::cout << "# 1D-spherical mapping." << std::endl;
        
        for (size_t j = 0; j != allvcom.size(); ++j) {
            vcom = allvcom[j];
            std::cout << "# " << vcom << std::endl;
            Tmt rngo(std::time(NULL));
            typedef boost::uniform_on_sphere<double, TVec > Tusph;
            boost::variate_generator<Tmt,Tusph > vgen( rngo,Tusph(3));
            
            size_t nbins = size_t(rmax/bin) ;
            
            TVec gofr(nbins,0);
            
            size_t c = 50000;
            for (size_t i = 0 ; i != nbins; ++i){
                double rt = i*bin;
                
                for (size_t k = 0; k != c; ++k ){
                    //gofr(i) += dens3d.get( rt * vgen() + vcom);
                    gofr(i) += dens3d.getper( rt * vgen() + vcom);
                }
                
                gofr(i) /= c;
                
                
            }
            
            for (size_t i = 0 ; i != nbins; ++i){
                std::cout << i*bin << " " << gofr(i) << std::endl;
            }
            
            std::cout << " " << std::endl;
        }
    } else if (maptype=="projxyz"){
        
        std::cout << "# 1D averaging along the x,y and z axes." << std::endl;
        
        // number of bins
        size_t nx =size_t(2*xmax/bin);
        size_t ny =size_t(2*ymax/bin);
        size_t nz =size_t(2*zmax/bin);
        
        
        // the grid must be centered on the solute
        //        double x0 = -1.0*xmax+vcom(0);
        //        double y0 = -1.0*ymax+vcom(1);
        //        double z0 = -1.0*zmax+vcom(2);
        
        TVec gridz(nz,0), gridy(ny,0), gridx(nx,0);
        
        dens3d.getProjXYZ(gridx,gridy,gridz,bulkdens);
        
        std::cout << "\n# z-axis" << std::endl;
        for (size_t i = 0; i != gridz.size();++i){
            std::cout << i << " " << gridz(i) << std::endl;
        }
        
        
    } else if (maptype=="excess"){
        
        dens3d.getIntExcess(bulkdens);
        
    } else if (maptype=="rhoel")
    {
        
        dens3d.getRhoEl();
        dens3d.writedxfile(filenameout);
        
      } else if (maptype=="rhoelreal") {

        dens3d.getRhoElReal();
        dens3d.writedxfile(filenameout);

    } else if (maptype=="blobs") {
        
        Tdens denshess(filenamehess,chemicalspecies);
        
        // condition moving fwd to having read the density file
        
        if (!denshess.getreadany()){
            std::cout << "[x] Input density file could not been read. Exiting..." << std::endl;
            return 1;
        }
        
        dens3d.blobs (denshess, threshold, bulkdens);
        
    }
    
    else if (maptype=="blobsper") {
        
        Tdens denshess(filenamehess,chemicalspecies);

        
        // condition moving fwd to having read the density file
        
        if (!denshess.getreadany()){
            std::cout << "[x] Input density file could not been read. Exiting..." << std::endl;
            return 1;
        }
        
        denshess.setBulkC(bulkdens);
        denshess.printDetails();
        dens3d.blobs_periodic(denshess, threshold);
        
      } else if (maptype=="cutresol"){


        std::cout << "Cutting density to keep a desired interval of resolutions." << std::endl;
        std::sort(resminmax.begin(),resminmax.end());
        std::cout << resminmax << std::endl;
        dens3d.cutResolution(resminmax[0],resminmax[1]);
        dens3d.writedxfile(filenameout);


      }




    else if (maptype=="string") {
        
        
        // where q, coms are stored
        std::vector< boost::math::quaternion <double> > vquaternion;
        std::vector< TVec > vcom;
        
        // (1) read the file and store com's and frames into vectors
        std::ifstream file(stringfile.c_str());
        std::string line;
        util::convto<double> tft;
        typedef std::vector<std::string> Tvs;
        Tvs stokens;
        
        
        size_t counter = 0;
        
        while( getline( file, line ) ) {
            
            stokens =  boost::copy_range<Tvs>(line |boost::adaptors::tokenized( std::string("[\\t\\s,]+"), -1 ) );
            
            if (counter == 0){
                if (stokens.size() == 3) {
                    
                    TVec v (3);
                    v(0) = tft(stokens[0]); v(1) = tft(stokens[1]); v(2) = tft(stokens[2]);
                    vcom.push_back(v);
                } else {
                    std::cout << " Something is worng with the worm." << std::endl;
                }
                counter = 1;
            } else {
                if (stokens.size() == 4) {
                    
                    boost::math::quaternion <double> q(tft(stokens[0]),tft(stokens[1]),tft(stokens[2]),tft(stokens[3]) ) ;
                    vquaternion.push_back(q);
                    
                } else {
                    std::cout << " Something is worng with the worm." << std::endl;
                }
                
                counter = 0;
            }
            
            
        }
        
        
        
        // (2) generate 3d points in the x,y plane, within the required distance from origin
        Tmt rngo(std::time(NULL));
        
        typedef boost::uniform_real <double> Tureal;
        
        boost::variate_generator<Tmt,Tureal > vgen( rngo, Tureal(0,3) );
        
        std::vector <  TVec  > xyplanevectors;
        for (size_t i = 0; i != 750; ++i) {
            
            TVec v(3,0);
            
            v(0) = vgen();
            v(1) = vgen();
            
            if (ublas::norm_2(v) <= 3.0)   xyplanevectors.push_back(v)  ;
        }
        
        std::cout <<"# Selected " << xyplanevectors.size() << " random vectors." << std::endl;
        
        double distance = 0.0;
        TVec origin(vcom[0]);
        std::vector <double> position;
        std::vector <double> density;
        
        for (size_t i = 0 ; i != ( vquaternion.size() - 1);++i)
        {
            
            size_t grain = 100;
            
            for (size_t j = 0; j != grain; ++j) {
                
                double t = j * (1.0f/double(grain));
                boost::math::quaternion<double> qintr;
                TVec comintr;
                rotations::slerp(vquaternion[i],vquaternion[i+1],t,qintr);
                comintr = (1.0f - t) * vcom[i] + t * vcom[i+1];
                
                
                distance += ublas::norm_2(comintr-origin);
                origin = comintr;
                position.push_back(distance);
                
                
                double dd = 0;
                
                for (size_t k = 0; k != xyplanevectors.size(); ++k) {
                    TVec v(3, 0);
                    v = rotations::rotatevector(xyplanevectors[k], qintr) + comintr;
                    dd += dens3d.get(v);
                }
                density.push_back(dd/xyplanevectors.size());
            }
            
            
        }
        
        
        for (size_t i = 0; i != density.size();++i)
        {
            
            std::cout << position[i] << " " << density[i] << std::endl;
            
        }
        
    } else if (maptype=="snake") {
        
        // read the snake file, series of 3D vectors
        // where q, coms are stored
        std::vector< TVec > vcom, vcominterpolated, vcominterpolatedder;
        
        // (1) read the file and store com's and frames into vectors
        std::ifstream file(stringfile.c_str());
        std::string line;
        util::convto<double> tft;
        typedef std::vector<std::string> Tvs;
        Tvs stokens;
        
        
        
        while( getline( file, line ) ) {
            
            stokens =  boost::copy_range<Tvs>(line |boost::adaptors::tokenized( std::string("[\\t\\s,]+"), -1 ) );
            
            if (stokens.size() == 3) {
                
                TVec v (3);
                v(0) = tft(stokens[0]); v(1) = tft(stokens[1]); v(2) = tft(stokens[2]);
                vcom.push_back(v);
            } else {
                std::cout << " Something is wrong with the worm. Expecting three numbers per line but got: " << line << std::endl;
            }
        }
        
        // do b-spline interpolation
        
        interpolation::CubicInterpolator< TVec > myInterpolator(100);
        
        
        myInterpolator.setCR();
        myInterpolator.setweightsf();
        
        
        
        myInterpolator.interpolate(vcom,vcominterpolated);
        
        
        
        std::cout << "# points along the string" << std::endl;
        
        double tmp = 0 ;
        
        for (size_t i = 0; i != vcominterpolated.size(); ++i) {
            
            if (i>0) tmp += ublas::norm_2(vcominterpolated[i-1]-vcominterpolated[i]);
            std::cout <<"# " << tmp << " "
            << vcominterpolated[i](0) << " "
            << vcominterpolated[i](1) << " "
            << vcominterpolated[i](2) << " "
            << std::endl ;
        }
        
        std::cout << " " << std::endl;
        
        myInterpolator.setweightsdf();
        
        myInterpolator.interpolate(vcom,vcominterpolatedder);
        
        
        // generates tangent to the curve by normalizing the spline derivative
        for (size_t i = 0; i != vcominterpolatedder.size(); ++i) {
            vcominterpolatedder[i] *= 1.0 / ublas::norm_2(vcominterpolatedder[i]);
        }
        
        
        
        // create a frame of reference, with points in the x,y plane
        
        // (2) generate 3d points in the x,y plane, within the required distance from origin
        
        double radius = rmax;
        
        
        Tmt rngo(std::time(NULL));
        
        typedef boost::uniform_real <double> Tureal;
        
        
        //boost::variate_generator<Tmt,Tureal > vgen( rngo, Tureal(0,2.0 * radius) );
        boost::variate_generator<Tmt,Tureal > vgen( rngo, Tureal(-radius,radius) );
        
        std::vector <  TVec  > xyplanevectors;
        for (size_t i = 0; i != 5000; ++i) {
            
            TVec v(3,0);
            
            //v(0) = vgen() - radius ;
            v(0) = vgen() ;
            //v(1) = vgen() - radius ;
            v(1) = vgen() ;
            
            if (ublas::norm_2(v) <= radius)   xyplanevectors.push_back(v)  ;
        }
        
        std::cout <<"# Selected " << xyplanevectors.size() << " random vectors." << std::endl;
        
        
        
        
        std::vector < TVec > xyplanevectorsquad;
        std::vector <double> wxyplanevectosrquad;
        
        double absc[10] = {0.0765265211334973,
            0.2277858511416451,
            0.3737060887154195,
            0.5108670019508271,
            0.6360536807265150,
            0.7463319064601508,
            0.8391169718222188,
            0.9122344282513259,
            0.9639719272779138,
            0.9931285991850949};
        
        double w[10] = {0.1527533871307258,
            0.1491729864726037,
            0.1420961093183820,
            0.1316886384491766,
            0.1181945319615184,
            0.1019301198172404,
            0.0832767415767048,
            0.0626720483341091,
            0.0406014298003869,
            0.0176140071391521};
        // r , theta
        
        
        // for each theta 0:360:resol
        //     for each r in abscisas.size
        //          xyplanevectorsquad = (r * cos(theta), r * sin(theta) )
        //          wxyplanevectosrquad = radius *  wgauss[i] * ( 2*pi*absc[i] )/resol
        
        
        for (size_t i = 0; i != 36; ++i){
            
            double theta = i * 2 * boost::math::constants::pi<double>()/36;
            
            for (size_t j = 0; j!= 10; ++j){
                
                TVec v(3,0);
                
                v(0) = absc[j] * radius * std::cos(theta);
                v(1) = absc[j] * radius * std::sin(theta);
                
                xyplanevectorsquad.push_back(v);
                
                wxyplanevectosrquad.push_back( radius * w[j] * 2.0 * boost::math::constants::pi<double>() * absc[j] * radius / 36 );
                
            }
            
        }
        
        
        double distance = 0.0;
        TVec origin(vcominterpolated[0]);
        std::vector <double> position;
        std::vector <double> density, densityquad;
        std::vector <double> densityvar;
        
        for (size_t i = 0; i != vcominterpolated.size(); ++i) {
            // find the rotation quaternion taking the predetermined points
            // a quaternion can be written easily in an axis-angle representation
            
            TVec axis(3,0);
            
            // projection of the tangent on the z-axis is tangent(2)
            // using cos(fi) = cos^2(fi/2) - sin^2(fi/2)
            double cosfihalf = std::sqrt(0.5* (  1 + vcominterpolatedder [i] (2)  ));
            double sinfihalf = std::sqrt(0.5* (  1 - vcominterpolatedder [i] (2)  ));
            
            // cross product of the tangent and z-axis
            axis (0)  = -vcominterpolatedder[i](1);
            axis (1)  =  vcominterpolatedder[i](0);
            axis (2)  =  0.0 ;
            // normalize axis
            // but making sure we don't normalize a null axis
            if (ublas::norm_2(axis)>1e-8) {
                axis *= (1.0 / ublas::norm_2(axis));
            }
            
            // multiply axis with sinfihalf
            axis *= sinfihalf;
            
            boost::math::quaternion<double> q(cosfihalf,axis(0),axis(1),axis(2));
            
            distance += ublas::norm_2(vcominterpolated[i]-origin);
            origin = vcominterpolated[i];
            position.push_back(distance);
            
            double dd = 0;
            
            
            std::vector<double> rotdens;
            
            for (size_t k = 0; k != xyplanevectors.size(); ++k) {
                
                TVec v(3, 0);
                
                v = rotations::rotatevector(xyplanevectors[k], q) + vcominterpolated[i] ;
                
                double temp  = dens3d.get(v);
                
                dd += temp;
                rotdens.push_back(temp);
            }
            
            
            density.push_back(dd/xyplanevectors.size());
            
            
            
            
            dd = 0;
            
            for (size_t k = 0 ; k != xyplanevectorsquad.size();++k){
                
                
                TVec v(3, 0);
                
                v = rotations::rotatevector(xyplanevectorsquad[k], q) + vcominterpolated[i] ;
                
                double temp  = dens3d.get(v);
                
                dd +=  (temp * wxyplanevectosrquad[k])  ;
                
                
                
                
            }
            
            
            densityquad.push_back(dd);
            
            
            
            
            
            
            
            if (dd > 1e-8) {
                
                TVec avevect(3, 0.0);
                
                
                for (size_t t = 0; t != xyplanevectors.size(); ++t) avevect = avevect + rotdens[t] * xyplanevectors[t];
                
                
                avevect /= dd;
                
                
                TVec varvect(3, 0.0);
                
                for (size_t t = 0; t != xyplanevectors.size(); ++t) {
                    
                    TVec tmp(3, 0);
                    
                    tmp = xyplanevectors[t] - avevect;
                    
                    tmp(0) = tmp(0) * tmp(0);
                    tmp(1) = tmp(1) * tmp(1);
                    tmp(2) = tmp(2) * tmp(2);
                    
                    varvect = varvect + rotdens[t] * tmp;
                }
                
                varvect /= dd;
                
                
                densityvar.push_back(0.5 * (std::sqrt(varvect(0)) + std::sqrt(varvect(1))));
                
            } else {
                //std::cout << ">>>>> " << i << " " << q <<  std::endl;
                densityvar.push_back(0.0);
            }
            
        }
        
        
        
        //for (size_t i = 0; i != density.size();++i)
        //{
        //    std::cout << position[i] << " " << density[i] * radius * radius * boost::math::constants::pi<double>() << std::endl;
        //}
        
        
        std::cout << " " << std::endl;
        
        
        
        //double multconst = bulkdens * 6.023 * 0.0001 *  radius * radius * boost::math::constants::pi<double>();
        
        for (size_t i = 0; i != density.size();++i)
        {
            std::cout << position[i] << " "
            //<< density[i] * multconst
            //<< " " << densityvar[i] * multconst
            << " " << densityquad[i] * 6.023 * 0.0001 * bulkdens
            << std::endl;
        }
        
        
        //std::cout << " " << std::endl;
        
        //for (size_t i = 0; i != densityvar.size();++i)
        //{
        //    std::cout << position[i] << " " << densityvar[i]  << std::endl;
        //}
        
        
        
    } // end snake approach
    
    
    else if (maptype=="worm") {
        
        // read the snake file, series of 3D vectors
        // where q, coms are stored
        std::vector< TVec > vcom, vcominterpolated, vcominterpolatedder;
        
        // (1) read the file and store com's and frames into vectors
        std::ifstream file(stringfile.c_str());
        std::string line;
        util::convto<double> tft;
        typedef std::vector<std::string> Tvs;
        Tvs stokens;
        
        
        
        while( getline( file, line ) ) {
            
            stokens =  boost::copy_range<Tvs>(line |boost::adaptors::tokenized( std::string("[\\t\\s,]+"), -1 ) );
            
            if (stokens.size() == 3) {
                
                TVec v (3);
                v(0) = tft(stokens[0]); v(1) = tft(stokens[1]); v(2) = tft(stokens[2]);
                vcom.push_back(v);
            } else {
                std::cout << " Something is wrong with the worm. Expecting three numbers per line but got: " << line << std::endl;
            }
        }
        
        
        
        
        // for ncycles
        
        
        std::vector <TVec> vcomp = vcom ;
        size_t ncontrolp =  32;
        
        for (size_t nc = 0; nc !=200; ++nc) {
            
            // do b-spline interpolation
            
            interpolation::CubicInterpolator<TVec> myInterpolator(100);
            
            //myInterpolator.setCR();
            //myInterpolator.setweightsf();
            
            
            myInterpolator.interpolate(vcomp, vcominterpolated);
            
            // compute the length of the curve segment
            
            
            double distance = 0.0;
            
            
            for (size_t i = 1; i != vcominterpolated.size(); ++i) {
                distance += ublas::norm_2(vcominterpolated[i] - vcominterpolated[i - 1]);
            }
            
            
            // chose 20 equidistant points along the curve
            double segment = distance/(ncontrolp - 1);
            
            std::vector<TVec> vcomequid;
            
            
            vcomequid.push_back(vcominterpolated.front());
            
            distance = 0.0;
            
            for (size_t i = 1; i != vcominterpolated.size() - 1; ++i) {
                
                distance += ublas::norm_2(vcominterpolated[i] - vcominterpolated[i - 1]);
                
                if ((distance > segment) && ( vcomequid.size() <= (ncontrolp - 2) )) {
                    
                    vcomequid.push_back(vcominterpolated[i - 1]);
                    distance = ublas::norm_2(vcominterpolated[i] - vcominterpolated[i - 1]);
                    
                }
                
            }
            
            vcomequid.push_back(vcominterpolated.back());
            
            
            // create a frame of reference (x,y,z points in a sphere)
            
            double radius = segment * 0.2;
            
            Tmt rngo(std::time(NULL));
            
            typedef boost::uniform_real<double> Tureal;
            
            boost::variate_generator<Tmt, Tureal> vgen(rngo, Tureal(0, radius * 2.0));
            
            std::vector<TVec> xyplanevectors;
            for (size_t i = 0; i != 750; ++i) {
                
                TVec v(3, 0);
                
                v(0) = vgen() - radius;
                v(1) = vgen() - radius;
                v(2) = vgen() - radius;
                
                if (ublas::norm_2(v) <= radius) xyplanevectors.push_back(v);
            }
            
            std::cout << "# Selected " << xyplanevectors.size() << " random vectors." << std::endl;
            
            
            std::vector<TVec> newvcomequid;
            
            for (size_t i = 0; i != vcomequid.size(); ++i) {
                
                
                TVec v = vcomequid[i];
                double dv = dens3d.get(v);
                
                std::cout << "# " << dv ;
                for (size_t k = 0; k != xyplanevectors.size(); ++k) {
                    
                    TVec vp = xyplanevectors[k] + vcomequid[i];
                    double dvp = dens3d.get(vp);
                    
                    if (dvp > dv) {
                        
                        dv = dvp;
                        
                        v = vp;
                    }
                }
                
                newvcomequid.push_back(v);
                std::cout << " vs " << dv << std::endl;
            }
            
            
            
            
            
            // add the first and the last points of the old list at the front and back of the new control points list
            
            
            newvcomequid.insert(newvcomequid.begin(), 1, vcom.front());
            newvcomequid.insert(newvcomequid.end(), 1, vcom.back());
            
            
            std::cout << " " << std::endl;
            
            for (size_t i = 0; i != newvcomequid.size(); ++i) {
                
                std::cout <<  newvcomequid[i] << std::endl;
            }
            
            
            
            vcomp = newvcomequid;
            
        }
        /*
         
         
         TVec origin(vcominterpolated[0]);
         std::vector <double> position;
         std::vector <double> density;
         std::vector <double> densityvar;
         
         for (size_t i = 0; i != vcominterpolated.size(); ++i) {
         double dd = 0;
         
         
         std::vector<double> rotdens;
         
         for (size_t k = 0; k != xyplanevectors.size(); ++k) {
         TVec v(3, 0);
         v = xyplanevectors[k] + vcominterpolated[i] ;
         dd += dens3d.get(v);
         rotdens.push_back(dens3d.get(v));
         }
         
         
         density.push_back(dd/xyplanevectors.size());
         
         
         TVec avevect(3,0.0);
         
         for (size_t t = 0; t != xyplanevectors.size() ; ++t) avevect = avevect + rotdens[t] * xyplanevectors[t];
         
         avevect /= dd;
         
         
         
         }
         
         
         */
        
    }
    
    
    //    else {
    //
    //
    //        std::cout << "Unknown mapping type (for the moment)." << std::endl;
    //
    //
    //    }
    
    
    
    std::cout << "# metatwist completed succesfully!"  << std::endl;
    
    
    return 0;
}
