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

#ifndef FFTWWRAP_HPP_
#define FFTWWRAP_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>


namespace ublas = boost::numeric::ublas;

#include "fftwmain.hpp"


namespace metafft {



template <typename T>
void FFT (std::vector< std::complex<T> > & data, std::vector< std::complex<T> > & ftdata, int forward);

template <typename T>
void FFT (ublas::vector< std::complex<T> > & data, ublas::vector< std::complex<T> > & ftdata, int forward);

template <typename T>
void FFT (ublas::vector<T> & data, ublas::vector< std::complex<T> > & ftdata, int forward);

template <typename T>
void FFT (boost::multi_array<std::complex<T>, 3> & data,boost::multi_array<std::complex<T>, 3> & ftdata, int forward );

template <typename T>
void FFT (boost::multi_array<T, 3> & data,boost::multi_array<std::complex<T>, 3> & ftdata, int forward );

template <typename T>
void FFT (ublas::vector<T> & data, ublas::vector<T> & ftdata, int forward);

template <typename T>
void FFT (ublas::vector<T> & data, ublas::vector< std::complex<T> > & ftdata, int forward){
	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.size() == ftdata.size() );


			if(forward == -1)  fftw::FFT1D_RC(&data.data()[0],&ftdata.data()[0],data.size());
			if(forward ==  1)  fftw::FFT1D_CR(&data.data()[0],&ftdata.data()[0],data.size());
}


template <typename T>
void FFT (std::vector< std::complex<T> > & data, std::vector< std::complex<T> > & ftdata, int forward){
	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.size() == ftdata.size() );

			if (forward == -1)  fftw::FFT1D_CC(&data.data()[0],&ftdata.data()[0],data.size(),forward);
			if (forward ==  1)  fftw::FFT1D_CC(&ftdata.data()[0],&data.data()[0],data.size(),forward);
}



template <typename T>
void FFT (ublas::vector< std::complex<T> > & data, ublas::vector< std::complex<T> > & ftdata, int forward){
	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.size() == ftdata.size() );

            if (forward == -1)  fftw::FFT1D_CC(&data.data()[0],&ftdata.data()[0],data.size(),forward);
            if (forward ==  1)  fftw::FFT1D_CC(&ftdata.data()[0],&data.data()[0],data.size(),forward);

}




template <typename T>
void FFT (boost::multi_array<std::complex<T>, 3> & data, boost::multi_array<std::complex<T>, 3> & ftdata, int forward ){
	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.shape()[0] == ftdata.shape()[0] );
	assert(data.shape()[1] == ftdata.shape()[1] );
	assert(data.shape()[2] == ftdata.shape()[2] );
			if (forward == -1) fftw::FFT3D_CC(&data.data()[0], &ftdata.data()[0], data.shape()[0], data.shape()[1], data.shape()[2], forward );
			if (forward ==  1) fftw::FFT3D_CC(&ftdata.data()[0], &data.data()[0], data.shape()[0], data.shape()[1], data.shape()[2], forward );
}



template <typename T>
void FFT (boost::multi_array<T, 3> & data,boost::multi_array<std::complex<T>, 3> & ftdata, int forward ){
// add
	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.shape()[0] == ftdata.shape()[0] );
	assert(data.shape()[1] == ftdata.shape()[1] );
	assert(data.shape()[2] >= int( float(ftdata.shape()[2]/2) + 1  ));
			if (forward == -1) fftw::FFT3D_RC(&data.data()[0],&ftdata.data()[0],data.shape()[0], data.shape()[1], data.shape()[2]) ;
			if (forward ==  1) fftw::FFT3D_CR(&data.data()[0],&ftdata.data()[0],data.shape()[0], data.shape()[1], data.shape()[2]) ;


}



template <typename T>
void FFT (ublas::vector<T> & data, ublas::vector<T> & ftdata, int forward){

	BOOST_STATIC_ASSERT((boost::is_same<T,double>::value) || (boost::is_same<T,float>::value));
	assert(data.size() == ftdata.size() );

	if(forward == -1)  fftw::FFT1D_DCTF(&data.data()[0],&ftdata.data()[0],data.size());
	if(forward ==  1)  fftw::FFT1D_DCTB(&data.data()[0],&ftdata.data()[0],data.size());
}




}
#endif /* FFTWWRAP_HPP_ */
