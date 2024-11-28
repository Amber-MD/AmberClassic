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

#ifndef FFTWMAIN_HPP_
#define FFTWMAIN_HPP_

#include <boost/static_assert.hpp>
#include <vector>
#include <cstring>
#include <complex>
#include <fftw3.h>

namespace metafft {
	namespace fftw {
//--------------------------------------------------
	template <typename T>
	void FFT1D_CC (std::complex<T> * data, std::complex<T> * ftdata, int size, int forward);

	template <>
	inline
	void FFT1D_CC (std::complex<double> * data, std::complex<double> * ftdata, int size, int forward);

	template <>
	inline
	void FFT1D_CC (std::complex<float> * data, std::complex<float> * ftdata, int size, int forward);

//--------------------------------------------------

	template <typename T>
	void FFT1D_RC (T * data, std::complex<T> * ftdata, int size);

	template <>
	inline
	void FFT1D_RC (double * data, std::complex<double> * ftdata, int size);

	template <>
	inline
	void FFT1D_RC (float * data, std::complex<float> * ftdata, int size);

//----------------------------------------------------

	template <typename T>
	void FFT1D_CR (T * data, std::complex<T> * ftdata, int size);

	template <>
	inline
	void FFT1D_CR (double * data, std::complex<double> * ftdata, int size);

	template <>
	inline
	void FFT1D_CR (float * data, std::complex<float> * ftdata, int size);


//----------------------------------------------------
	template <typename T>
	void FFT1D_DCTF(T * data, T * ftdata, int size);

	template <>
	inline
	void FFT1D_DCTF(float * data, float * ftdata, int size);

	template <>
	inline
	void FFT1D_DCTF(double * data, double * ftdata, int size);


	template <typename T>
	void FFT1D_DCTB(T * data, T * ftdata, int size);

	template <>
	inline
	void FFT1D_DCTB(float * data, float * ftdata, int size);

	template <>
	inline
	void FFT1D_DCTB(double * data, double * ftdata, int size);




//----------------------------------------------------
	template <typename T>
	void FFT3D_CC (std::complex<T> * data, std::complex<T> * ftdata, int size1, int size2, int size3, int forward);

	template <>
	inline
	void FFT3D_CC (std::complex<double> * data, std::complex<double> * ftdata, int size1, int size2, int size3, int forward);

	template <>
	inline
	void FFT3D_CC (std::complex<float> * data, std::complex<float> * ftdata, int size1, int size2, int size3, int forward);

//----------------------------------------------------
	template <typename T>
	void FFT3D_CR (T * data, std::complex<T> * ftdata, int size1, int size2, int size3);

	template <>
	inline
	void FFT3D_CR (double * data, std::complex<double> * ftdata, int size1, int size2, int size3);

	template <>
	inline
	void FFT3D_CR (float * data, std::complex<float> * ftdata, int size1, int size2, int size3);
//----------------------------------------------------
	template <typename T>
	void FFT3D_RC (T * data, std::complex<T> * ftdata, int size1, int size2, int size3);

	template <>
	inline
	void FFT3D_RC (double * data, std::complex<double> * ftdata, int size1, int size2, int size3);

	template <>
	inline
	void FFT3D_RC (float * data, std::complex<float> * ftdata, int size1, int size2, int size3);

//---------------FFT1d_CC----------------


//	template <typename T>
//	void FFT1D_CC (std::complex<T> * data, std::complex<T> * ftdata, int size, int forward){
//      BOOST_STATIC_ASSERT(false);
//	}

	template <>
	inline
	void FFT1D_CC (std::complex<double> * data, std::complex<double> * ftdata, int size, int forward){

		fftw_complex *pdata, *pftdata;
		pdata =   reinterpret_cast<fftw_complex*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);
		fftw_plan fwd;
		fwd  = fftw_plan_dft_1d(size, pdata, pftdata, forward, FFTW_ESTIMATE);
		fftw_execute(fwd);
		fftw_destroy_plan(fwd);

	}


	template <>
	inline
	void FFT1D_CC (std::complex<float> * data, std::complex<float> * ftdata, int size, int forward){

		fftwf_complex *pdata, *pftdata;
		pdata =   reinterpret_cast<fftwf_complex*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);
		fftwf_plan fwd;
		fwd  = fftwf_plan_dft_1d(size, pdata, pftdata, forward, FFTW_ESTIMATE);
		fftwf_execute(fwd);
		fftwf_destroy_plan(fwd);

	}

//-----------------FFT1D_RC---------------------------------


//	template <typename T>
//	void FFT1D_RC (T * data, std::complex<T> * ftdata, int size){
//		 BOOST_STATIC_ASSERT(false);
//	}


	template <>
	inline
	void FFT1D_RC (double * data, std::complex<double> * ftdata, int size){

		double *pdata;
		fftw_complex *pftdata;
		pdata =   reinterpret_cast<double*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);

		fftw_plan fwd;
		fwd  = fftw_plan_dft_r2c_1d(size, pdata, pftdata,FFTW_ESTIMATE);
	  	fftw_execute(fwd);
	  	fftw_destroy_plan(fwd);
	}

	template <>
	inline
	void FFT1D_RC (float * data, std::complex<float> * ftdata, int size){

		float *pdata;
		fftwf_complex *pftdata;
		pdata =   reinterpret_cast<float*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);

		fftwf_plan fwd;
		fwd  = fftwf_plan_dft_r2c_1d(size, pdata, pftdata,FFTW_ESTIMATE);
	  	fftwf_execute(fwd);
	  	fftwf_destroy_plan(fwd);
	}



//-----------------FFT1D_CR-----------------------------------------

//	template <typename T>
//	void FFT1D_CR (T * data, std::complex<T> * ftdata, int size){
//		BOOST_STATIC_ASSERT(false);
//
//	}

	template <>
	inline
	void FFT1D_CR (double * data, std::complex<double> * ftdata, int size){

		double *pdata;
		fftw_complex *pftdata;
		pdata =   reinterpret_cast<double*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);

		fftw_plan bkwd;
	  	bkwd = fftw_plan_dft_c2r_1d(size, pftdata, pdata,FFTW_ESTIMATE);
	 	fftw_execute(bkwd);
	  	fftw_destroy_plan(bkwd);

	}

	template <>
	inline
	void FFT1D_CR (float * data, std::complex<float> * ftdata, int size){

		float *pdata;
		fftwf_complex *pftdata;
		pdata =   reinterpret_cast<float*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);

		fftwf_plan bkwd;
	  	bkwd = fftwf_plan_dft_c2r_1d(size, pftdata, pdata,FFTW_ESTIMATE);
	 	fftwf_execute(bkwd);
	  	fftwf_destroy_plan(bkwd);

	}



	template <>
	inline
	void FFT1D_DCTF(double * data, double * ftdata, int size){

		double *pdata;
		double *pfdata;
		pdata = reinterpret_cast<double*>(data);
		pfdata = reinterpret_cast<double*>(ftdata);

		fftw_plan fwd;
		fwd = fftw_plan_r2r_1d(size, pdata, pfdata,FFTW_REDFT10, FFTW_ESTIMATE);
	 	fftw_execute(fwd);
	  	fftw_destroy_plan(fwd);

	}

	template <>
	inline
	void FFT1D_DCTF(float * data, float * ftdata, int size){

		float *pdata;
		float *pfdata;
		pdata = reinterpret_cast<float*>(data);
		pfdata = reinterpret_cast<float*>(ftdata);

		fftwf_plan fwd;
		fwd = fftwf_plan_r2r_1d(size, pdata, pfdata,FFTW_REDFT10, FFTW_ESTIMATE);
	 	fftwf_execute(fwd);
	  	fftwf_destroy_plan(fwd);

	}



	template <>
	inline
	void FFT1D_DCTB(double * data, double * ftdata, int size){

		double *pdata;
		double *pfdata;
		pdata = reinterpret_cast<double*>(data);
		pfdata = reinterpret_cast<double*>(ftdata);

		fftw_plan fwd;
		fwd = fftw_plan_r2r_1d(size, pfdata, pdata,FFTW_REDFT01, FFTW_ESTIMATE);
	 	fftw_execute(fwd);
	  	fftw_destroy_plan(fwd);

	}

	template <>
	inline
	void FFT1D_DCTB(float * data, float * ftdata, int size){

		float *pdata;
		float *pfdata;
		pdata = reinterpret_cast<float*>(data);
		pfdata = reinterpret_cast<float*>(ftdata);

		fftwf_plan fwd;
		fwd = fftwf_plan_r2r_1d(size, pfdata, pdata,FFTW_REDFT01, FFTW_ESTIMATE);
	 	fftwf_execute(fwd);
	  	fftwf_destroy_plan(fwd);

	}





//----------------------------FFT3D_CC----------------------------------------
//	template <typename T>
//	void FFT3D_CC(std::complex<T> * data, std::complex<T> * ftdata, int size1, int size2, int size3, int forward){
//		BOOST_STATIC_ASSERT(false);
//
//	}

	template <>
	inline
	void FFT3D_CC(std::complex<double> * data, std::complex<double> * ftdata, int size1, int size2, int size3, int forward){


		fftw_complex *pdata, *pftdata;
		pdata =   reinterpret_cast<fftw_complex*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);
		fftw_plan fwd;
		fwd  = fftw_plan_dft_3d(size1,size2,size3, pdata, pftdata, forward, FFTW_ESTIMATE);
		fftw_execute(fwd);
		fftw_destroy_plan(fwd);

	}

	template <>
	inline
	void FFT3D_CC(std::complex<float> * data, std::complex<float> * ftdata, int size1, int size2, int size3, int forward){


		fftwf_complex *pdata, *pftdata;
		pdata =   reinterpret_cast<fftwf_complex*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);
		fftwf_plan fwd;
		fwd  = fftwf_plan_dft_3d(size1,size2,size3, pdata, pftdata, forward, FFTW_ESTIMATE);
		fftwf_execute(fwd);
		fftwf_destroy_plan(fwd);

	}
//-------------------------------------------------------------------------------------------------------------
//	template <typename T>
//	void FFT3D_CR (T * data, std::complex<T> * ftdata, int size1, int size2, int size3){
//	// add
//		BOOST_STATIC_ASSERT(false);
//	}


	template <>
	inline
	void FFT3D_CR (double * data, std::complex<double> * ftdata, int size1, int size2, int size3){
	// add

		double *pdata;
		fftw_complex *pftdata;
		pdata =   reinterpret_cast<double*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);

		fftw_plan bkwd;
		bkwd  = fftw_plan_dft_c2r_3d(size1, size2, size3, pftdata, pdata,FFTW_ESTIMATE);
	  	fftw_execute(bkwd);
	  	fftw_destroy_plan(bkwd);
	}

	template <>
	inline
	void FFT3D_CR (float * data, std::complex<float> * ftdata, int size1, int size2, int size3){
	// add

		float *pdata;
		fftwf_complex *pftdata;
		pdata =   reinterpret_cast<float*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);

		fftwf_plan bkwd;
		bkwd  = fftwf_plan_dft_c2r_3d(size1, size2, size3, pftdata, pdata,FFTW_ESTIMATE);
	  	fftwf_execute(bkwd);
	  	fftwf_destroy_plan(bkwd);
	}


//------------------------------------------------------------------------------------------------------
//	template <typename T>
//	void FFT3D_RC (T * data, std::complex<T> * ftdata, int size1, int size2, int size3){
//	// add
//		BOOST_STATIC_ASSERT(false);
//	}



	template <>
	inline
	void FFT3D_RC (double * data, std::complex<double> * ftdata, int size1, int size2, int size3){
	// add

		double *pdata;
		fftw_complex *pftdata;
		pdata =   reinterpret_cast<double*>(data);
		pftdata = reinterpret_cast<fftw_complex*>(ftdata);

		fftw_plan fwd;
		fwd  = fftw_plan_dft_r2c_3d(size1, size2, size3, pdata, pftdata,FFTW_ESTIMATE);
	  	fftw_execute(fwd);
	  	fftw_destroy_plan(fwd);


	}

	template <>
	inline
	void FFT3D_RC (float * data, std::complex<float> * ftdata, int size1, int size2, int size3){
	// add

		float *pdata;
		fftwf_complex *pftdata;
		pdata =   reinterpret_cast<float*>(data);
		pftdata = reinterpret_cast<fftwf_complex*>(ftdata);

		fftwf_plan fwd;
		fwd  = fftwf_plan_dft_r2c_3d(size1, size2, size3, pdata, pftdata,FFTW_ESTIMATE);
	  	fftwf_execute(fwd);
	  	fftwf_destroy_plan(fwd);


	}




	}
}

#endif /* FFTWMAIN_HPP_ */
