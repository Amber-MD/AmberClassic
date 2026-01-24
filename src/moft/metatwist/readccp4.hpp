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


#ifndef READCCP4_HPP_
#include "boosthead.hpp"
#include <fstream>
#include <iostream>


namespace util {

  namespace files {

    template <typename T>
    bool writeccp4(const std::string & filename,
                   const ublas::vector <T> &orig,
                   const ublas::vector <T> &a, const ublas::vector <T> &b, const ublas::vector <T> &c,
                   const boost::multi_array<T, 3> &data3d){

      std::ofstream myFile (filename, std::ios::out | std::ios::binary);

      union headerbyte { unsigned char c[224]; int i[56]; float f[56]; };

      headerbyte headerccp;

      for (int i = 0; i != 56; ++i) headerccp.f[i] = (float)0.0;


      T radpi = 180.0/boost::math::constants::pi<T>();

      /* 1      NC (slowest)*/ headerccp.i[ 0] = data3d.shape()[2] ;
      /* 2      NR    v     */ headerccp.i[ 1] = data3d.shape()[1] ;
      /* 3      NS (fastest)*/ headerccp.i[ 2] = data3d.shape()[0] ;
      /* 4      MODE        */ headerccp.i[ 3] = 2  ;
      /* 5      NCSTART     */ headerccp.i[ 4] = 0  ;
      /* 6      NRSTART     */ headerccp.i[ 5] = 0  ;
      /* 7      NSSTART     */ headerccp.i[ 6] = 0  ;
      /* 8      NX          */ headerccp.i[ 7] = data3d.shape()[0] ;
      /* 9      NY          */ headerccp.i[ 8] = data3d.shape()[1] ;
      /*10      NZ          */ headerccp.i[ 9] = data3d.shape()[2] ;
      /*11      X length    */ headerccp.f[10] = ublas::norm_2(a) ;
      /*12      Y length    */ headerccp.f[11] = ublas::norm_2(b) ;
      /*13      Z length    */ headerccp.f[12] = ublas::norm_2(c) ;
      /*14      Alpha       */ headerccp.f[13] = radpi*std::acos(ublas::inner_prod(b,c)/(ublas::norm_2(b)*ublas::norm_2(c)));
      /*15      Beta        */ headerccp.f[14] = radpi*std::acos(ublas::inner_prod(a,c)/(ublas::norm_2(a)*ublas::norm_2(c)));
      /*16      Gamma       */ headerccp.f[15] = radpi*std::acos(ublas::inner_prod(b,a)/(ublas::norm_2(b)*ublas::norm_2(a)));
      /*17      MAPC        */ headerccp.i[16] = 3;
      /*18      MAPR        */ headerccp.i[17] = 2;
      /*19      MAPS        */ headerccp.i[18] = 1;
      /*20      AMIN        */ headerccp.f[19] = *std::min_element(data3d.data(),data3d.data()+data3d.num_elements());
      /*21      AMAX        */ headerccp.f[20] = *std::max_element(data3d.data(),data3d.data()+data3d.num_elements());
      /*22      AMEAN       */ headerccp.f[21] = std::accumulate(data3d.data(),data3d.data()+data3d.num_elements(),0.0)/data3d.num_elements();
      /*23      ISPG        */ headerccp.i[22] = 1;
      /*24      NSYMBT      */ headerccp.i[23] = 0;
      /*25      LSKFLG      */ headerccp.i[24] = 0;
      /*26-34   SKWMAT(S11) */ headerccp.f[25] = 0.0;
      /*26-34   SKWMAT(S12) */ headerccp.f[26] = 0.0;
      /*26-34   SKWMAT(S13) */ headerccp.f[27] = 0.0;
      /*26-34   SKWMAT(S21) */ headerccp.f[28] = 0.0;
      /*26-34   SKWMAT(S22) */ headerccp.f[29] = 0.0;
      /*26-34   SKWMAT(S23) */ headerccp.f[30] = 0.0;
      /*26-34   SKWMAT(S31) */ headerccp.f[31] = 0.0;
      /*26-34   SKWMAT(S32) */ headerccp.f[32] = 0.0;
      /*26-34   SKWMAT(S33) */ headerccp.f[33] = 0.0;
      /*35-37   SKWTRN(T1)  */ headerccp.f[34] = 0.0;
      /*35-37   SKWTRN(T2)  */ headerccp.f[35] = 0.0;
      /*35-37   SKWTRN(T3)  */ headerccp.f[36] = 0.0;
      /*38-52   future use  */
      /*50-52  origin       */ headerccp.f[49] = orig(0) ; headerccp.f[50] = orig(1) ; headerccp.f[51] = orig(2);
      /*53      MAP         */ headerccp.c[4*52]   = 'M'; headerccp.c[4*52+1] = 'A';
      headerccp.c[4*52+2] = 'P'; headerccp.c[4*52+3] = ' ';

      /*54      Machine Stmp*/ headerccp.c[4*53]   = 'D' ; headerccp.c[4*53+1] = 'A' ;
      headerccp.c[4*53+2] = '\0'; headerccp.c[4*53+3] = '\0' ;


      // summarize header info
      std::cout<<"# |   Summarizing the CPP4 file features for:" << std::endl;
      std::cout<<"# |   "<< filename << std::endl;
      std::cout<<"# |   (based upon description found at https://www.ccpem.ac.uk/mrc_format/mrc2014.php)." << std::endl;
      std::cout<<"# |   1      NC (slowest)"<< headerccp.i[ 0]             << std::endl;
      std::cout<<"# |   2      NR    v     "<< headerccp.i[ 1]             << std::endl;
      std::cout<<"# |   3      NS (fastest)"<< headerccp.i[ 2]             << std::endl;
      std::cout<<"# |   4      MODE        "<< headerccp.i[ 3]             << std::endl;
      std::cout<<"# |   5      NCSTART     "<< headerccp.i[ 4]             << std::endl;
      std::cout<<"# |   6      NRSTART     "<< headerccp.i[ 5]             << std::endl;
      std::cout<<"# |   7      NSSTART     "<< headerccp.i[ 6]             << std::endl;
      std::cout<<"# |   8      NX          "<< headerccp.i[ 7]             << std::endl;
      std::cout<<"# |   9      NY          "<< headerccp.i[ 8]             << std::endl;
      std::cout<<"# |  10      NZ          "<< headerccp.i[ 9]             << std::endl;
      std::cout<<"# |  11      X length    "<< headerccp.f[10]             << std::endl;
      std::cout<<"# |  12      Y length    "<< headerccp.f[11]             << std::endl;
      std::cout<<"# |  13      Z length    "<< headerccp.f[12]             << std::endl;
      std::cout<<"# |  14      Alpha       "<< headerccp.f[13]             << std::endl;
      std::cout<<"# |  15      Beta        "<< headerccp.f[14]             << std::endl;
      std::cout<<"# |  16      Gamma       "<< headerccp.f[15]             << std::endl;
      std::cout<<"# |  17      MAPC        "<< headerccp.i[16]             << std::endl;
      std::cout<<"# |  18      MAPR        "<< headerccp.i[17]             << std::endl;
      std::cout<<"# |  19      MAPS        "<< headerccp.i[18]             << std::endl;
      std::cout<<"# |  20      AMIN        "<< headerccp.f[19]             << std::endl;
      std::cout<<"# |  21      AMAX        "<< headerccp.f[20]             << std::endl;
      std::cout<<"# |  22      AMEAN       "<< headerccp.f[21]             << std::endl;
      std::cout<<"# |  23      ISPG        "<< headerccp.i[22]             << std::endl;
      std::cout<<"# |  24      NSYMBT      "<< headerccp.i[23]             << std::endl;
      std::cout<<"# |  25      LSKFLG      "<< headerccp.i[24]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S11) "<< headerccp.f[25]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S12) "<< headerccp.f[26]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S13) "<< headerccp.f[27]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S21) "<< headerccp.f[28]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S22) "<< headerccp.f[29]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S23) "<< headerccp.f[30]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S31) "<< headerccp.f[31]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S32) "<< headerccp.f[32]             << std::endl;
//      std::cout<<"# |  26-34   SKWMAT(S33) "<< headerccp.f[33]             << std::endl;
//      std::cout<<"# |  35-37   SKWTRN(T1)  "<< headerccp.f[34]             << std::endl;
//      std::cout<<"# |  35-37   SKWTRN(T2)  "<< headerccp.f[35]             << std::endl;
//      std::cout<<"# |  35-37   SKWTRN(T3)  "<< headerccp.f[36]             << std::endl;
      std::cout<<"# |  50      ORIGIN      "<< headerccp.f[49]              << std::endl;
      std::cout<<"# |  51      ORIGIN      "<< headerccp.f[50]              << std::endl;
      std::cout<<"# |  52      ORIGIN      "<< headerccp.f[51]              << std::endl;
      //std::cout<<" 38-52   future use  "/*<< headerccp.i[ ] */         << std::endl;
      std::cout<<"# |  53      MAP         "<< headerccp.c[4*52]<<headerccp.c[4*52+1]<<headerccp.c[4*52+2]<<headerccp.c[4*52+3]<< std::endl;
      std::cout<<"# |  54      MACHST      "<< boost::format("0x%02x0x%02x0x%02x0x%02x")
                 % (int) headerccp.c[4*53]
                 % (int) headerccp.c[4*53+1]
          % (int) headerccp.c[4*53+2]
          % (int) headerccp.c[4*53+3]
          << std::endl;
      std::cout<<"# |  55      ARMS        "<< headerccp.f[55]             << std::endl;

      // write the header
      myFile.write(reinterpret_cast<char*>(&headerccp.c[0]), 224*sizeof(unsigned char));

      // write the empty "for future use" section.
      float labels[200];

      for (size_t i = 0 ; i != 200; ++i) labels[i]=(float)0.0;

      myFile.write(reinterpret_cast<char*>(&labels[0]), sizeof(labels) );

      size_t s = sizeof (float);
      for (size_t i = 0 ; i != data3d.num_elements();++i){
          float t = data3d.data()[i];
          myFile.write(reinterpret_cast<char*>(&t),s);
        }
      myFile.close();


      return true;
    }



    template <typename T>
    bool readccp4(const std::string &filename,
                  ublas::vector <T> &orig,
                  ublas::matrix <T> &rot,
                  boost::multi_array<T, 3> &data3d){

      /*
             (from http://www.ccp4.ac.uk/html/maplib.html)
             The overall layout of the file is as follows:
             
             File header (256 longwords)
             Symmetry information
             Map, stored as a 3-dimensional array
             
             The header is organised as 56 words followed by space for ten 80 character text labels as follows:
             1       NC              # of Columns    (fastest changing in map)
             2       NR              # of Rows
             3       NS              # of Sections   (slowest changing in map)
             4       MODE            Data type
                        0 = envelope stored as signed bytes (from
                        -128 lowest to 127 highest)
                        1 = Image     stored as Integer*2
                        2 = Image     stored as Reals (Mode 2 is the normal mode used in CCP4.)
                        3 = Transform stored as Complex Integer*2
                        4 = Transform stored as Complex Reals
                        5 == 0
             5       NCSTART         Number of first COLUMN  in map
             6       NRSTART         Number of first ROW     in map
             7       NSSTART         Number of first SECTION in map
             8       NX              Number of intervals along X
             9       NY              Number of intervals along Y
             10      NZ              Number of intervals along Z
             11      X length        Cell Dimensions (Angstroms)
             12      Y length                     "
             13      Z length                     "
             14      Alpha           Cell Angles     (Degrees)
             15      Beta                         "
             16      Gamma                        "
             17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
             18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
             19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
             20      AMIN            Minimum density value
             21      AMAX            Maximum density value
             22      AMEAN           Mean    density value    (Average)
             23      ISPG            Space group number
             24      NSYMBT          Number of bytes used for storing symmetry operators
             25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
             26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
             LSKFLG .ne. 0.
             35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
             Skew transformation is from standard orthogonal
             coordinate frame (as used for atoms) to orthogonal
             map frame, as
             
             Xo(map) = S * (Xo(atoms) - t)
             
             38      future use       (some of these are used by the MSUBSX routines
             .          "              in MAPBRICK, MAPCONT and FRODO)
             .          "   (all set to zero by default)
             .          "
             52          "
             
             53    MAP            Character string 'MAP ' to identify file type
             54    MACHST        Machine stamp indicating the machine type
             which wrote file
             55      ARMS            Rms deviation of map from mean density
             56      NLABL           Number of labels being used
             57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)
             */


      std::ifstream myFile (filename, std::ios::in | std::ios::binary);



      if (myFile) {

          //std::cout << ">> Reading from: " << filename << std::endl;

          union headerbyte { unsigned char c[224]; int i[56]; float f[56]; };
          headerbyte headerccp;

          // read header
          myFile.read((char*)&headerccp.i, 56*sizeof(int));
          // read labels
          char labels[801];
          myFile.read(labels, 800*(sizeof(char)));



          // summarize header info
          std::cout<<"# |   Summarizing the CPP4 file features for:" << std::endl;
          std::cout<<"# |   "<< filename << std::endl;
          std::cout<<"# |   (based upon description found at https://www.ccpem.ac.uk/mrc_format/mrc2014.php)." << std::endl;
          std::cout<<"# |   1      NC (slowest)"<< headerccp.i[ 0]             << std::endl;
          std::cout<<"# |   2      NR    v     "<< headerccp.i[ 1]             << std::endl;
          std::cout<<"# |   3      NS (fastest)"<< headerccp.i[ 2]             << std::endl;
          std::cout<<"# |   4      MODE        "<< headerccp.i[ 3]             << std::endl;
          std::cout<<"# |   5      NCSTART     "<< headerccp.i[ 4]             << std::endl;
          std::cout<<"# |   6      NRSTART     "<< headerccp.i[ 5]             << std::endl;
          std::cout<<"# |   7      NSSTART     "<< headerccp.i[ 6]             << std::endl;
          std::cout<<"# |   8      NX          "<< headerccp.i[ 7]             << std::endl;
          std::cout<<"# |   9      NY          "<< headerccp.i[ 8]             << std::endl;
          std::cout<<"# |  10      NZ          "<< headerccp.i[ 9]             << std::endl;
          std::cout<<"# |  11      X length    "<< headerccp.f[10]             << std::endl;
          std::cout<<"# |  12      Y length    "<< headerccp.f[11]             << std::endl;
          std::cout<<"# |  13      Z length    "<< headerccp.f[12]             << std::endl;
          std::cout<<"# |  14      Alpha       "<< headerccp.f[13]             << std::endl;
          std::cout<<"# |  15      Beta        "<< headerccp.f[14]             << std::endl;
          std::cout<<"# |  16      Gamma       "<< headerccp.f[15]             << std::endl;
          std::cout<<"# |  17      MAPC        "<< headerccp.i[16]             << std::endl;
          std::cout<<"# |  18      MAPR        "<< headerccp.i[17]             << std::endl;
          std::cout<<"# |  19      MAPS        "<< headerccp.i[18]             << std::endl;
          std::cout<<"# |  20      AMIN        "<< headerccp.f[19]             << std::endl;
          std::cout<<"# |  21      AMAX        "<< headerccp.f[20]             << std::endl;
          std::cout<<"# |  22      AMEAN       "<< headerccp.f[21]             << std::endl;
          std::cout<<"# |  23      ISPG        "<< headerccp.i[22]             << std::endl;
          std::cout<<"# |  24      NSYMBT      "<< headerccp.i[23]             << std::endl;
          std::cout<<"# |  25      LSKFLG      "<< headerccp.i[24]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S11) "<< headerccp.f[25]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S12) "<< headerccp.f[26]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S13) "<< headerccp.f[27]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S21) "<< headerccp.f[28]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S22) "<< headerccp.f[29]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S23) "<< headerccp.f[30]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S31) "<< headerccp.f[31]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S32) "<< headerccp.f[32]             << std::endl;
//          std::cout<<"# |  26-34   SKWMAT(S33) "<< headerccp.f[33]             << std::endl;
//          std::cout<<"# |  35-37   SKWTRN(T1)  "<< headerccp.f[34]             << std::endl;
//          std::cout<<"# |  35-37   SKWTRN(T2)  "<< headerccp.f[35]             << std::endl;
//          std::cout<<"# |  35-37   SKWTRN(T3)  "<< headerccp.f[36]             << std::endl;
          std::cout<<"# |  50      ORIGIN      "<< headerccp.f[49]             << std::endl;
          std::cout<<"# |  51      ORIGIN      "<< headerccp.f[50]             << std::endl;
          std::cout<<"# |  52      ORIGIN      "<< headerccp.f[51]             << std::endl;
          //std::cout<<" 38-52   future use  "/*<< headerccp.i[ ] */         << std::endl;
          std::cout<<"# |  53      MAP         "<< headerccp.c[4*52]<<headerccp.c[4*52+1]<<headerccp.c[4*52+2]<<headerccp.c[4*52+3]<< std::endl;
          std::cout<<"# |  54      MACHST      "<< boost::format("0x%02x0x%02x0x%02x0x%02x")
                                                                % (int) headerccp.c[4*53]
                                                                % (int) headerccp.c[4*53+1]
                                                                % (int) headerccp.c[4*53+2]
                                                                % (int) headerccp.c[4*53+3]
                                                                << std::endl;
//          std::cout<<"# |  55      ARMS        "<< headerccp.f[55]             << std::endl;
          //std::cout<<" 56      NLABL       "<< headerccp.i[56]             << std::endl;
          //std::cout<<" 57-256  LABEL(20,10)\n"<< Labels         << std::endl;


          // mapping between cartesian and grid axes.
          ublas::vector<size_t> gr2ax(3), ax2gr(3);

          gr2ax(0) = headerccp.i[16]-1 ;
          gr2ax(1) = headerccp.i[17]-1 ;
          gr2ax(2) = headerccp.i[18]-1 ;

          ax2gr(gr2ax(0)) = 0;
          ax2gr(gr2ax(1)) = 1;
          ax2gr(gr2ax(2)) = 2;

          // figure out the cell vectors
          // a || x ; (a,b) gamma ; (a,c) beta ; (b,c) alpha
          T pirad = boost::math::constants::pi<T>()/180.0;
          T alpha = pirad *  headerccp.f[13];
          T beta  = pirad *  headerccp.f[14];
          T gamma = pirad *  headerccp.f[15];

          T ca = std::cos(alpha);
          T cb = std::cos(beta);
          T cg = std::cos(gamma);




          // make sure rot has the right size
          rot.resize(3,3);

          // from cosines to unit cell vectors
          rot(0,0) = 1.0;
          rot(1,0) = 0.0;
          rot(2,0) = 0.0;

          rot(0,1) = cg ;
          rot(1,1) = std::sqrt ( 1.0 - cg*cg  );
          rot(2,1) = 0.0;

          rot(0,2) = cb ;
          rot(1,2) = (ca-cg*cb) / rot(1,1);
          rot(2,2) = std::sqrt ( 1.0 - rot(0,2)*rot(0,2) - rot(1,2)*rot(1,2));

          rot(0,0) *= ( headerccp.f[10] / headerccp.i[7] );
          rot(1,0) *= ( headerccp.f[10] / headerccp.i[7] );
          rot(2,0) *= ( headerccp.f[10] / headerccp.i[7] );
          rot(0,1) *= ( headerccp.f[11] / headerccp.i[8] );
          rot(1,1) *= ( headerccp.f[11] / headerccp.i[8] );
          rot(2,1) *= ( headerccp.f[11] / headerccp.i[8] );
          rot(0,2) *= ( headerccp.f[12] / headerccp.i[9] );
          rot(1,2) *= ( headerccp.f[12] / headerccp.i[9] );
          rot(2,2) *= ( headerccp.f[12] / headerccp.i[9] );


          orig.resize(3,0.0);

          if ( headerccp.f[49]==0.0 and headerccp.f[50]==0.0 and headerccp.f[51]==0.0 ){
            orig(0) = headerccp.i[ 4 + gr2ax(0) ];
            orig(1) = headerccp.i[ 4 + gr2ax(1) ];
            orig(2) = headerccp.i[ 4 + gr2ax(2) ];
            orig =  ublas::prod(rot,orig);
            std::cout << "# Origin NxNyNz: " << orig << std::endl;
          } else {
            orig(0) = headerccp.f[49];
            orig(1) = headerccp.f[50];
            orig(2) = headerccp.f[51];
            // std::cout << "# Origin MRC [v2000]: " << orig << std::endl;
          };

          // Populating the tmp data array:
          typename boost::multi_array<T, 3>::size_type ordering[] = {ax2gr(0),ax2gr(1),ax2gr(2)};
          bool ascending[] = {true,true,true};
          boost::general_storage_order<3> so(ordering,ascending);


          boost::multi_array<T, 3> tmp3d(boost::extents[ headerccp.i[ax2gr(0)] ]
                                                       [ headerccp.i[ax2gr(1)] ]
                                                       [ headerccp.i[ax2gr(2)] ]
                                                       , so
                                                        );



          for (size_t i = 0 ; i != tmp3d.num_elements();++i){
                 float t;
                 myFile.read(reinterpret_cast<char*>(&t),sizeof(float));
                 tmp3d.data()[i] = t;
            }


          if (data3d.num_elements()==0){


              data3d.resize(boost::extents[ headerccp.i[ax2gr(0)] ]
                                          [ headerccp.i[ax2gr(1)] ]
                                          [ headerccp.i[ax2gr(2)] ]);

              for (size_t i = 0; i != data3d.num_elements(); ++i) data3d.data()[i] = 0.0;


            } else if (data3d.num_elements() != headerccp.i[2] * headerccp.i[1] * headerccp.i[0]) {

              std::exit(0);

            } else {

              std::cout << "# Accumulating data." << std::endl;

            }



          size_t idim = data3d.shape()[0];
          size_t jdim = data3d.shape()[1];
          size_t kdim = data3d.shape()[2];


              for (size_t i = 0; i != idim ; ++i) {
                  for (size_t j = 0; j != jdim ; ++j) {
                      for (size_t k = 0; k != kdim ; ++k) {
                          data3d[i][j][k] += tmp3d[i][j][k];
                        }
                    }
                }

          return true;

        } else {



          std::cout << ">>  Could not read file:" << filename << std::endl;
          return false;
        }


      myFile.close();


    }

  }
}

#endif /* READCCP4_HPP_ */

