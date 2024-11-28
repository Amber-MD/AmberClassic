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


#ifndef CONVCORR_HPP_
#define CONVCORR_HPP_

#include "fftwwrap.hpp"

namespace metafft{

//--------Convolutions------

template <typename T>
void convolution( ublas::vector< std::complex<T> > & g,  ublas::vector<std::complex<T> > & h, ublas::vector< std::complex<T> > & corr);

template <typename T>
void convolution( ublas::vector< T > & g,  ublas::vector< T > & h, ublas::vector< T > & corr);

template <typename T>
void convolution( boost::multi_array<std::complex<T>, 3>  & g,  boost::multi_array<std::complex<T>, 3> & h, boost::multi_array<std::complex<T>, 3> & corr);



template <typename T>
void convolution( boost::multi_array<std::complex<T>, 3>  & g,  boost::multi_array<std::complex<T>, 3> & h);



//---------Correlation Facility--------------//

template <typename T>
void correlation( ublas::vector< T > & g,  ublas::vector< T > & h, ublas::vector< T > & corr, size_t lag);

template <typename T>
void correlation( std::vector< ublas::vector< T > > & g, std::vector< ublas::vector< T > > & h, std::vector < ublas::vector< T > > & corr, size_t lag);

template <typename T>
void correlation( std::vector< ublas::vector< T > > & g, std::vector< ublas::vector< T > > & h, std::vector < ublas::vector< T > > & corr);


//template <typename T>
//void correlation( std::vector< ublas::vector< std::complex<T> > > & g, std::vector< ublas::vector< std::complex<T> > > & h, std::vector < ublas::vector< std::complex<T> > > & corr);

template <typename T>
void correlation( ublas::vector< std::complex<T> > & g,  ublas::vector<std::complex<T> > & h, ublas::vector< std::complex<T> > & corr);

template <typename T>
void correlation( ublas::vector< T > & g,  ublas::vector< T > & h, ublas::vector< T > & corr);




//------Implementation-------------------------//


template <typename T>
void correlation( std::vector< ublas::vector< T > > & g,
		          std::vector< ublas::vector< T > > & h,
		          std::vector < ublas::vector< T > > & corr){

	corr.resize(g.size());
	for (size_t i = 0; i != g.size(); ++i)	correlation(g[i],h[i], corr[i]);

}


template <typename T>
void correlation( std::vector< ublas::vector< T > > & g,
		          std::vector< ublas::vector< T > > & h,
		          std::vector < ublas::vector< T > > & corr, size_t lag){

	corr.resize(g.size());
	for (size_t i = 0; i < g.size(); ++i)	{
		correlation(g[i],h[i],corr[i], lag);
	}

}


template <typename T>
void correlation( ublas::vector< T > & g,  ublas::vector< T > & h, ublas::vector< T > & corr, size_t lag){

	bool preserve = true;

	size_t s = g.size() + lag;
	g.resize(s, preserve);
	h.resize(s, preserve);


	for (size_t i = g.size()-lag; i != g.size(); ++i){
		g(i) = 0;
		h(i) = 0;
	}


	correlation(g,h,corr);
}










template <typename T>
void correlation( ublas::vector< std::complex<T> > & g,  ublas::vector<std::complex<T> > & h, ublas::vector< std::complex<T> > & corr){

		corr.resize(g.size());
		corr *= 0.0;

		ublas::vector< std::complex<T> > gft(g.size()), hft(h.size()), corrft(corr.size());
		corrft *= 0.0;
		hft *= 0.0;
		gft *= 0.0;



		FFT(h,hft,-1);
		FFT(g,gft,-1);
		for (size_t i = 0; i < gft.size(); ++i) corrft(i) = gft(i)*std::conj(hft(i));
	    FFT(corr,corrft,1);
	    corr /= g.size();

}


template <typename T>
void correlation( ublas::vector< T > & g,  ublas::vector<T > & h, ublas::vector< T > & corr){

		corr.resize(g.size());
		corr *= 0.0;

		ublas::vector< std::complex<T> > gft(g.size()), hft(h.size()), corrft(corr.size());
		corrft *= 0.0;
		hft *= 0.0;
		gft *= 0.0;

		FFT(h,hft,-1);
		FFT(g,gft,-1);
		for (size_t i = 0; i < gft.size(); ++i) corrft(i) = gft(i)*std::conj(hft(i));
	    FFT(corr,corrft,1);
	    corr /= g.size();


}








template <typename T>
void convolution( boost::multi_array<T, 3>  & g,  boost::multi_array<T, 3> & h, boost::multi_array<T, 3> & conv) {

    typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
    typedef boost::multi_array<T, 3> Tarray;
    typedef typename Tarray::index index;

    // make sure g and h have the same size
    assert(g.shape()[2] == h.shape()[2] );
    assert(g.shape()[1] == h.shape()[1] );
    assert(g.shape()[0] == h.shape()[0] );


    index xdim = g.shape()[0];
    index ydim = g.shape()[1];
    index zdim = g.shape()[2];




    boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
    boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, (zdim/2)+1}};

    // resize conv
    conv.resize(shape);

    Tarray_c gfft (shapeft), hfft(shapeft);

    metafft::FFT(g,gfft,-1);
    metafft::FFT(h,hfft,-1);

    for (size_t i = 0; i != hfft.num_elements();++i){
        gfft.data()[i] *= hfft.data()[i];
    }

    metafft::FFT(conv,gfft,1);

    T norm = T(conv.num_elements());

    for (size_t i = 0; i != conv.num_elements();++i){
        conv.data()[i] /= norm;
    }


}


template <typename T>
void convolution( boost::multi_array<T, 3>  & g,  boost::multi_array<T, 3> & h) {

    typedef boost::multi_array<std::complex<T>, 3> Tarray_c;
    typedef boost::multi_array<T, 3> Tarray;
    typedef typename Tarray::index index;

    // make sure g and h have the same size
    assert(g.shape()[2] == h.shape()[2] );
    assert(g.shape()[1] == h.shape()[1] );
    assert(g.shape()[0] == h.shape()[0] );


    index xdim = g.shape()[0];
    index ydim = g.shape()[1];
    index zdim = g.shape()[2];

//    boost::array<typename Tarray::index, 3> shape = {{ xdim, ydim, zdim}};
    boost::array<typename Tarray::index, 3> shapeft = {{ xdim, ydim, (zdim/2)+1}};

    Tarray_c gfft (shapeft), hfft(shapeft);

    metafft::FFT(g,gfft,-1);
    metafft::FFT(h,hfft,-1);

    for (size_t i = 0; i != hfft.num_elements();++i){
        gfft.data()[i] *= hfft.data()[i];
    }

    metafft::FFT(g,gfft,1);

    T norm = T(g.num_elements());

    for (size_t i = 0; i != g.num_elements();++i){
        g.data()[i] /= norm;
    }

}







}
#endif /* CONVCORR_HPP_ */
