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









#ifndef MOFT_INTERPOLATION_HPP
#define MOFT_INTERPOLATION_HPP


#include <utilities/boosthead.hpp>
#include <utilities/utilities.hpp>



namespace interpolation{



    template <typename T>
    class CubicInterpolator {


    public:

        CubicInterpolator(size_t grain): grain(grain)
        {

            this->matrix.resize(4,4);
            this->weights.resize(this->grain);

            this->setBS();
            this->setweightsf();

            std::cout <<"# Initializing cubic interpolator." << std::endl;

        };

        //  set weights for computing function
        void setweightsf(){


            ublas::vector<double> t(this->grain,0.0);

            double tspacing = 1.0/double(this->grain);

            for (size_t i = 0; i != t.size(); ++i){
                t(i) = i * tspacing;
            }

            for (size_t j = 0; j != t.size(); ++j){
                double tt = t(j);
                ublas::vector< double > uvec(4,1.0);
                uvec(0) = tt * tt * tt;
                uvec(1) = tt * tt;
                uvec(2) = tt;
                //uvec(3)=1.0;
                this->weights[j] = ublas::prec_prod(uvec,this->matrix);
            }

        }

        //  set weights for computing first derivative
        void setweightsdf(){


            ublas::vector<double> t(this->grain,0.0);

            double tspacing = 1.0/double(this->grain);

            for (size_t i = 0; i != t.size(); ++i){
                t(i) = i * tspacing;
            }

            for (size_t j = 0; j != t.size(); ++j){
                double tt = t(j);
                ublas::vector< double > uvec(4,0.0);
                uvec(0) = 3.0 * tt * tt;
                uvec(1) = 2.0 * tt;
                uvec(2) = 1.0 ;
                //uvec(3)=1.0;
                this->weights[j] = ublas::prec_prod(uvec,this->matrix);
            }

        }



        // interpolate based on the aforementioned weights
        void interpolate( const std::vector< T > & input, std::vector< T> & output ){

            output.resize(0);

            for (size_t i = 1; i != input.size()-2; ++i) {
                for (size_t j = 0; j !=  this->weights.size(); ++j ){

                    output.push_back(

                           this->weights[j](0) * input[i-1]
                         + this->weights[j](1) * input[i]
                         + this->weights[j](2) * input[i+1]
                         + this->weights[j](3) * input[i+2]

                    );


                }

            }

        };

        //  set Catmull-Rom Interpolation
        void setCR(){

            this->matrix(0,0)=-1.0; this->matrix(0,1)= 3.0; this->matrix(0,2)=-3.0; this->matrix(0,3)= 1.0;
            this->matrix(1,0)= 2.0; this->matrix(1,1)=-5.0; this->matrix(1,2)= 4.0; this->matrix(1,3)=-1.0;
            this->matrix(2,0)=-1.0; this->matrix(2,1)= 0.0; this->matrix(2,2)= 1.0; this->matrix(2,3)= 0.0;
            this->matrix(3,0)= 0.0; this->matrix(3,1)= 2.0; this->matrix(3,2)= 0.0; this->matrix(3,3)= 0.0;

            this->matrix *= 0.5;

        };
        //  set cubic B-spline interpolation
        void setBS(){

            this->matrix(0,0)=-1.0; this->matrix(0,1)= 3.0; this->matrix(0,2)=-3.0; this->matrix(0,3)= 1.0;
            this->matrix(1,0)= 3.0; this->matrix(1,1)=-6.0; this->matrix(1,2)= 3.0; this->matrix(1,3)= 0.0;
            this->matrix(2,0)=-3.0; this->matrix(2,1)= 0.0; this->matrix(2,2)= 3.0; this->matrix(2,3)= 0.0;
            this->matrix(3,0)= 1.0; this->matrix(3,1)= 4.0; this->matrix(3,2)= 1.0; this->matrix(3,3)= 0.0;
            this->matrix *= 1.0/6.0;

        };

    private:


        size_t grain;
        ublas::matrix<double> matrix;
        std::vector <  ublas::vector <double> > weights;


    };
}
#endif //MOFT_INTERPOLATION_HPP
