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


#ifndef TRANSROT_HPP_
#define TRANSROT_HPP_

#include "utilities.hpp"
#include<algorithm>


template <typename T>
class TransMat{

	typedef ublas::matrix<T> Tmat;
public:
	TransMat(){
		mat.resize(3,4,0);
		for (size_t i = 0; i < 3; ++i)	mat(i,i)=1.0;
	};

const	Tmat & get(){return mat;}

	void set(Tmat & xmat) {

		assert(xmat.size1()==3);
		assert(xmat.size2()==4);

		mat.data() = xmat.data();
	}

private:

	Tmat mat;

};

/**
 *
 * @param coor
 * @param trans
 */
template <typename T>
void TransMatTransf(ublas::matrix<T> & coor, TransMat<T> & trans){

	ublas::matrix<T>  coortmp(coor);
	coortmp.resize(coortmp.size1(),4,true);
	for (size_t i = 0; i != coortmp.size1(); ++i)  coortmp(i,3) = 1.0 ;
	for (size_t i = 0; i != coor.size1(); ++i) {
	 ublas::row(coor,i) = ublas::prec_prod(trans.get(),  ublas::row(coortmp,i)  );
	}
}

/**
 *
 * @param coors
 * @param refcoors
 * @param tr
 */
template <typename T>
void rmsd(ublas::matrix<T> & coors, ublas::matrix<T> & refcoors, TransMat<T> & tr ){

	assert(coors.size2()==3);
	assert(refcoors.size2()==3);
	assert(refcoors.size1()==coors.size1());


	ublas::vector<T> centers(3,0);
	ublas::vector<T> refcenters(3,0);


	for (size_t i = 0; i < coors.size2() ; ++i){
		centers(i) = ublas::sum(ublas::column(coors,i) )/coors.size1() ;
		refcenters(i) = ublas::sum(ublas::column(refcoors,i) )/refcoors.size1() ;
	}

	for (size_t i = 0; i != coors.size1(); ++i){
		ublas::row(coors,i) -= centers;
		ublas::row(refcoors,i) -= refcenters;
	}

	ublas::matrix<T> R(3,3,0);

	for (size_t i = 0; i < 3 ; ++i){
		for (size_t j = 0; j < 3 ; ++j){
			const ublas::matrix_column< ublas::matrix<T> > t1(refcoors,j),t2(coors,i);
			R(i,j)= ublas::prec_inner_prod(t1,t2);
		}
	}


	ublas::matrix<T,ublas::column_major> S(4,4,0);

	  S(0, 0) = R(0, 0) + R(1, 1) + R(2, 2);
	  S(1, 0) = R(1, 2) - R(2, 1);
	  S(2, 0) = R(2, 0) - R(0, 2);
	  S(3, 0) = R(0, 1) - R(1, 0);

	  S(0, 1) = S(1, 0);
	  S(1, 1) = R(0, 0) - R(1, 1) - R(2, 2);
	  S(2, 1) = R(0, 1) + R(1, 0);
	  S(3, 1) = R(0, 2) + R(2, 0);

	  S(0, 2) = S(2, 0);
	  S(1, 2) = S(2, 1);
	  S(2, 2) =-R(0, 0) + R(1, 1) - R(2, 2);
	  S(3, 2) = R(1, 2) + R(2, 1);

	  S(0, 3) = S(3, 0);
	  S(1, 3) = S(3, 1);
	  S(2, 3) = S(3, 2);
	  S(3, 3) =-R(0, 0) - R(1, 1) + R(2, 2);

		ublas::vector<T> Sv;

		diagfirsteigen(S,Sv);

		ublas::matrix<T> U(3,4,0);

		rotation_matrix(U,Sv);
		//U(3,3) = 1.0;
		U(0,3) = refcenters(0)-centers(0);
		U(1,3) = refcenters(1)-centers(1);
		U(2,3) = refcenters(2)-centers(2);

		tr.set(U);
}

/**
 *
 * @param A
 * @param eigvc
 */

//template<typename T>
//void diagfirsteigen(ublas::matrix<T,ublas::column_major > & A, ublas::vector<T> & eigvc){
//
//	assert(A.size1()==A.size2());
//	size_t n = A.size1();
//    ublas::vector<std::complex<T> > values(n);
//    ublas::matrix<std::complex<T>, ublas::column_major> Vectors_left(n,n);
//    ublas::matrix<std::complex<T>, ublas::column_major> Vectors_right(n,n);
//
//    ulapack::geev(A, values, &Vectors_left, &Vectors_right, ulapack::optimal_workspace());
//
//    size_t max = 0;
//    for (size_t j = 0 ; j != values.size(); ++j) if( values(j).real() > values(max).real() ) max = j;
//
//    eigvc.resize(n);
//    for(size_t i = 0; i < Vectors_right.size1(); ++i) eigvc(i) = Vectors_right(i,max).real();
//
//}

/**
 *
 * @param U
 * @param q
 */
template<typename T>
void rotation_matrix(ublas::matrix<T> & U, ublas::vector<T> & q){

	T q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33;

	  q0 = q(0);
	  q1 = q(1);
	  q2 = q(2);
	  q3 = q(3);

	  b0 = 2.0 * q0;
	  b1 = 2.0 * q1;
	  b2 = 2.0 * q2;
	  b3 = 2.0 * q3;

	  q00 = b0 * q0-1.0;
	  q01 = b0 * q1;
	  q02 = b0 * q2;
	  q03 = b0 * q3;

	  q11 = b1 * q1;
	  q12 = b1 * q2;
	  q13 = b1 * q3;

	  q22 = b2 * q2;
	  q23 = b2 * q3;

	  q33 = b3 * q3;

	  U(0,0) = q00 + q11;
	  U(0,1) = q12 - q03;
	  U(0,2) = q13 + q02;

	  U(1,0) = q12 + q03;
	  U(1,1) = q00 + q22;
	  U(1,2) = q23 - q01;

	  U(2,0) = q13 - q02;
	  U(2,1) = q23 + q01;
	  U(2,2) = q00 + q33;

}

#endif /* TRANSROT_HPP_ */
