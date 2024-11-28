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



#ifndef LUTILITIES_HPP_
#define LUTILITIES_HPP_


namespace util {

template<class T>
bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
	typedef ublas::permutation_matrix<std::size_t> pmatrix;
	inverse.resize(input.size1(),input.size2());
	// create a working copy of the input
	ublas::matrix<T> A(input);
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
	// perform LU-factorization
	int res = ublas::lu_factorize(A,pm);
	if( res != 0 ) {
		return false;
	} else {
		// create identity matrix of "inverse"
		inverse.assign(ublas::identity_matrix<T>(A.size1()));
		// backsubstitute to get the inverse
		ublas::lu_substitute(A, pm, inverse);
		return true;
	}
}


//
//size_t getIndex(const size_t & g1,const size_t & g2,const size_t & g3,const size_t & n1, const size_t & n2, const size_t & n3){
//
//	return g1*n2*n3 + g2*n3 + g3;
//
//}



template<typename T>
float blend103(const  T & r,    const  T &  s,    const  T &  t,   const  T & x000, const T & x001,
	       const  T & x010, const  T &  x011, const  T & x100, const  T & x101,
	       const  T & x110, const  T &  x111)
{

        return  x000
                + r         * ( - x000 + x100 )
                +     s     * ( - x000        + x010 )
                +         t * ( - x000               + x001 )
                + r * s     * (   x000 - x100 - x010                      + x110 )
                + r     * t * (   x000 - x100        - x001        + x101 )
                +     s * t * (   x000        - x010 - x001 + x011 )
                + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 );

}
}



#endif /* LUTILITIES_HPP_ */
