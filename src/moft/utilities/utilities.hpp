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


#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>


#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH
#include <boost/foreach.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <utilities/statistics.hpp>
//#include <utilities/lapack.hpp>

namespace ublas = boost::numeric::ublas;
//namespace po = boost::program_options;
namespace bl = boost::lambda;



namespace std {
template <typename T>
std::istream& operator>>(std::istream& in, std::pair<T,T>& p) {
    in >> p.first >> p.second;
    return in;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::pair<T,T>& p) {
    out << "("<< p.first << ","<< p.second << ")\n";
    return out;
}


}


template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
   return os;
}



template<typename T>
std::ostream& operator<<(std::ostream &out, const std::complex<T> &rhs){
  out << (T)rhs.real() << " " << (T)rhs.imag();
  return out;
}


template <typename T>
bool operator==(const ublas::matrix<T> & a, const ublas::matrix<T> & b){

	if ( (a.size1()!=b.size1()) || (a.size2()!=b.size2())  ) return false;

	for (size_t i = 0; i != a.data().size(); ++i){
		if(a.data()[i] != b.data()[i]  ) return false;
	}
	return true;
}



template <typename T>
bool operator==(const ublas::vector<T> & a, const ublas::vector<T> & b){

	if ( (a.size()!=b.size())   ) {
		return false;
	}

	for (size_t i = 0; i != a.data().size(); ++i){
		if(a.data()[i] != b.data()[i]  ) {
			return false;
		}
	}
    return true;
}








namespace util {
/**
 *
 * @param pdbfilename
 * @param names
 * @param numbers
 */
template <typename T>
void readPDB(const std::string & pdbfilename, std::vector< std::string > & names, std::vector<T> & numbers, std::vector<std::string> & residues);

/**
 *
 * @param pdbfilename
 * @param number
 * @param name
 * @param resname
 * @param chain
 * @param resnum
 * @param occ
 * @param beta
 */
template <typename T>
void readPDB3(const std::string & pdbfilename
			   , std::vector<size_t> & number
			   , std::vector<std::string> & name
			   , std::vector<std::string> & resname
			   , std::vector<std::string> & chain
			   , std::vector<size_t> & resnum
			   , std::vector<T> & occ
			   , std::vector<T> & beta
//		           , std::vector<std::string> & symbol
//		           , std::vector<T> & charge
);


/**
 *
 * @param X
 * @param Sigma
 */
template <typename T>
void bartletterror(const ublas::vector<T> & X,
                                 ublas::vector<T> & Sigma);
/**
 *
 * @param v
 */
template <typename T>
void invert_square (ublas::vector<T> & v);

/**
 *
 * @param mat
 * @param lambda
 */
template <typename T>
void scale_diag(ublas::matrix<T> & mat, T lambda);

/**
 *
 * @param mat
 * @param filename
 */
template <typename T>
void read(ublas::matrix<T> & mat, std::string filename);

/**
 *
 * @param vec
 * @param filename
 */
template <typename T>
void read(ublas::vector<T> & vec, std::string filename);



/**
 *
 * @param mat
 * @param filename
 */
template <typename T>
void write(const ublas::matrix<T> & mat, const std::string & filename);

/**
 *
 * @param vec
 * @param filename
 */
template <typename T>
void write(const ublas::vector<T> & vec, const std::string & filename);
/**
 *
 * @param vec
 * @param filename
 */
template <typename T>
void write(const std::vector< ublas::vector<T> > & vec, const std::string & filename);

template <typename T>
void writetimeser (const std::vector< ublas::vector<T> > & vec, const float & step, const std::string & filename );


template <typename T>
size_t round(T x)
{
        //std::cout <<" round " << x << std::endl;
   return size_t(x > 0.0 ? x + 0.5 : x - 0.5);
        //return x > 0.0 ? x + 0.5 : x - 0.5;
}

/**
 * converts anyting to a anything else
 */

template <typename T>
struct convto

{
	typedef T result_type;
	template <typename TI>
	T operator()(const TI & s ) const {
		return boost::lexical_cast<T> (s);
	};

};




template <typename T>
void readPDB3(const std::string & pdbfilename
			   , std::vector<size_t> & number
			   , std::vector<std::string> & name
			   , std::vector<std::string> & resname
			   , std::vector<std::string> & chain
			   , std::vector<size_t> & resnum
			   , std::vector<T> & occ
			   , std::vector<T> & beta
//		           , std::vector<std::string> & symbol
//		           , std::vector<T> & charge
)


{
//	COLUMNS        DATA  TYPE    FIELD        DEFINITION
//	-------------------------------------------------------------------------------------
//	 1 -  6        Record name   "ATOM  "
//	 7 - 11        Integer       serial       Atom  serial number.
//	13 - 16        Atom          name         Atom name.
//	17             Character     altLoc       Alternate location indicator.
//	18 - 20        Residue name  resName      Residue name.
//	22             Character     chainID      Chain identifier.
//	23 - 26        Integer       resSeq       Residue sequence number.
//	27             AChar         iCode        Code for insertion of residues.
//	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
//	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
//	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
//	55 - 60        Real(6.2)     occupancy    Occupancy.
//	61 - 66        Real(6.2)     tempFactor   Temperature  factor.
//	77 - 78        LString(2)    element      Element symbol, right-justified.
//	79 - 80        LString(2)    charge       Charge  on the atom.

using boost::trim;
using boost::lexical_cast;

	std::ifstream xfile(pdbfilename.c_str());
	if(! xfile) {
		std::cout << pdbfilename << " is invalid" << std::endl;
		throw;
	}
	std::string line;
	while( getline( xfile, line ) ) {

		if(line.substr(0,6)=="ATOM  "){
			std::string temp;

			temp = line.substr(6,5);
			trim(temp);
			number.push_back(lexical_cast<size_t>(temp));

			temp = line.substr(12,4);
			trim(temp);
			name.push_back(temp);

			temp = line.substr(17,4);
			trim(temp);
			resname.push_back(temp);

			temp = line.substr(21,1);
			trim(temp);
			chain.push_back(temp);

			temp = line.substr(22,4);
			trim(temp);
			resnum.push_back(lexical_cast<size_t>(temp));

			temp = line.substr(54,6);
			trim(temp);
			occ.push_back(lexical_cast<T>(temp));

			temp = line.substr(60,6);
			trim(temp);
			beta.push_back(lexical_cast<T>(temp));

		}
	}








}


template <typename T>
void readPDB(const std::string & pdbfilename, std::vector< std::string > & names, std::vector<T> & numbers, std::vector<std::string> & residues){


	typedef boost::tokenizer< boost::char_separator<char> > tokenizer ;
	typedef tokenizer::iterator tokenizerit;

	std::ifstream xfile(pdbfilename.c_str());
	if(! xfile) {
		std::cout << pdbfilename << " is invalid" << std::endl;
		throw;
	}
	std::string line;
	boost::char_separator<char> sep(" ");

	std::cout << ">> Reading pdb data from "<< pdbfilename << " file ." << std::endl;

	while( getline( xfile, line ) ) {
		tokenizer tok(line,sep);
		tokenizerit beg=tok.begin();
		if(boost::lexical_cast<std::string>(*beg++)=="ATOM"){
			numbers.push_back(boost::lexical_cast<T>(*beg++));
			names.push_back(boost::lexical_cast<std::string>(*beg++));
			residues.push_back(boost::lexical_cast<std::string>(*beg++));
		};
	}

	xfile.close();

	if(( names.size()!= numbers.size() ) || ( names.size()!= residues.size() ) || (residues.size()!= numbers.size())) {
		throw;
	} else {
		std::cout << ">> Read "<< pdbfilename << " containing " << numbers.size() << " records." << std::endl;
	}
}


template <typename T>
void read(ublas::matrix<T> & mat, const std::string filename){



	typedef std::vector<std::string> Tvs;
	Tvs stokens;

	std::ifstream file(filename.c_str());
	std::string line;

	size_t linecounter = 0, dim = 0, len = 0;

	std::cout << ">> Reading n-dimensional (n>=2) from "<< filename << " file ." << std::endl;

	while( getline( file, line ) ) {

		if (linecounter == 0){
            stokens =  boost::copy_range<Tvs>(line |boost::adaptors::tokenized( std::string("[\\t\\s,]+"), -1 ) );
			len = boost::lexical_cast<int>(stokens[1]);
			dim = boost::lexical_cast<int>(stokens[2]);
			mat.resize(len,dim);
			std::cout << ">> Data size is " << len << " x "<< dim <<"." << std::endl;
		} else if (linecounter <= len) {
			//std::cout << "-------> "<< line << std::endl;
			boost::algorithm::trim(line);
			//std::cout << "-------> "<< line << std::endl;
            boost::copy(line | boost::adaptors::tokenized( std::string("[\\t\\s,]+"), -1 )
			     | boost::adaptors::transformed( convto<double>() ),
				     mat.data().begin()+dim*(linecounter-1)
				 );
		}
		++linecounter;
	}


	file.close();


}




template <typename T>
void read(ublas::vector<T> & vec,  std::string filename){

	typedef std::vector<std::string> Tvs;
	Tvs stokens;

	std::ifstream file(filename.c_str());
	std::string line;


	size_t linecounter = 0, len = 0;


	std::cout << ">> Reading from "<< filename << " file ." << std::endl;

	while( getline( file, line ) ) {
		if (linecounter == 0){
			stokens =  boost::copy_range<Tvs>(line|boost::adaptors::tokenized( "[\\t\\s,]+", -1 ) );
			vec.resize( boost::lexical_cast<int>(stokens[1]) );
			std::cout << ">> Data size is " << len << "." << std::endl;
		} else if (linecounter <= len) {
			stokens =  boost::copy_range<Tvs>(line|boost::adaptors::tokenized( "[\\t\\s,]+", -1 ) );
			vec(linecounter-1) = boost::lexical_cast<T>(stokens[0]);
		}
		++linecounter;
	}

	file.close();
}







template <typename T>
void write(const ublas::matrix<T> & mat, const std::string & filename){

	std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;
	file << "# "<<mat.size1()<< " " << mat.size2() << "\n";

	for (size_t i = 0; i < mat.size1(); ++i){

		for (size_t j = 0; j < mat.size2(); ++j){

			file << mat(i,j) << " ";


		}

		file << "\n";
	}

	file.close();




}



template <typename T>
void writepath(const ublas::matrix<T> & mat, const std::string & path, const std::string & filename){

	boost::filesystem::path p("./"+path);

	if ( !boost::filesystem::is_directory(p) ) {
		boost::filesystem::create_directories(p);
		std::cout << ">> Created directory: " << path << std::endl;
	} else {
		std::cout << ">> Directory allready exists: " << path << std::endl;

	}
	boost::filesystem::path q(path+"/"+filename);
	boost::filesystem::ofstream file;
	file.open(q);


	//std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;
	file << "# "<<mat.size1()<< " " << mat.size2() << "\n";

	for (size_t i = 0; i < mat.size1(); ++i){

		for (size_t j = 0; j < mat.size2(); ++j){

			file << mat(i,j) << " ";


		}

		file << "\n";
	}

	file.close();




}





template <typename T>
void writepath(const ublas::vector<T> & vec, const std::string & path, const std::string & filename){


	boost::filesystem::path p("./"+path);

	if ( !boost::filesystem::is_directory(p) ) {
		boost::filesystem::create_directories(p);
		std::cout << ">> Created directory: " << path << std::endl;
	} else {
		std::cout << ">> Directory allready exists: " << path << std::endl;

	}
	boost::filesystem::path q(path+"/"+filename);
	boost::filesystem::ofstream file;
	file.open(q);



	//std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;
	file <<"# " <<vec.size()<< "\n";

	for (size_t i = 0; i < vec.size(); ++i){

			file << vec(i) << "\n";

	}

	file.close();

}








template <typename T>
void write(const ublas::vector<T> & vec, const std::string & filename){

	std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;
	file <<"# " <<vec.size()<< "\n";

	for (size_t i = 0; i < vec.size(); ++i){

			file << vec(i) << "\n";

	}

	file.close();

}





template <typename T>
void write(const std::vector< ublas::vector<T> > & vec, const std::string & filename){

	std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;

	file <<"# " << vec.size() << " " << vec[0].size() << "\n";

	for (size_t i = 0; i < vec.size(); ++i){
	   for (size_t j = 0; j != vec[0].size(); ++j){


	   // std::cout << vec[i](j) << " ";
			file << vec[i](j) << " ";

		}
		 file << "\n";
	 //std::cout << std::endl;
	}

	file.close();
}



template <typename T>
void writetimeser (const std::vector< ublas::vector<T> > & vec, const float & step,const std::string & filename ) {

	std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;

	file <<"# " << vec[0].size() << " " << vec.size() +1 << "\n";

	for (size_t j = 0; j != vec[0].size(); ++j){
		file << j*step << " ";
		for (size_t i = 0; i < vec.size(); ++i){
			file << vec[i](j) << " ";
		}
		 file << "\n";
	}

	file.close();


}


template <typename T>
void writetimeserpath (const std::vector< ublas::vector<T> > & vec, const float & step,const std::string & path,  const std::string & filename ) {

	boost::filesystem::path p("./"+path);

	if ( !boost::filesystem::is_directory(p) ) {
		boost::filesystem::create_directories(p);
		std::cout << ">> Created directory: " << path << std::endl;
	} else {
		std::cout << ">> Directory allready exists: " << path << std::endl;

	}
	boost::filesystem::path q(path+"/"+filename);
	boost::filesystem::ofstream file;
	file.open(q);



	//std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;

	file <<"# " << vec[0].size() << " " << vec.size() +1 << "\n";

	for (size_t j = 0; j != vec[0].size(); ++j){
		file << j*step << " ";
		for (size_t i = 0; i < vec.size(); ++i){
			file << vec[i](j) << " ";
		}
		 file << "\n";
	}

	file.close();


}


template <typename T>
void writetimeserpath (const  ublas::vector<T>  & vec, const float & step,const std::string & path,  const std::string & filename ) {

	boost::filesystem::path p("./"+path);

	if ( !boost::filesystem::is_directory(p) ) {
		boost::filesystem::create_directories(p);
		std::cout << ">> Created directory: " << path << std::endl;
	} else {
		std::cout << ">> Directory allready exists: " << path << std::endl;

	}
	boost::filesystem::path q(path+"/"+filename);
	boost::filesystem::ofstream file;
	file.open(q);



	//std::ofstream file(filename.c_str());

	std::cout << "<< Writing to "<< filename << " file ." << std::endl;

	file <<"# " << vec.size() << "\n";

	for (size_t j = 0; j != vec.size(); ++j){
		file << j*step << " "<< vec[j];
		file << "\n";
	}

	file.close();


}




}
#endif /* UTILITIES_HPP_ */
