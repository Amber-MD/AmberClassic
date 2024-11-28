#ifndef LAPACK_HPP
#define LAPACK_HPP



#include <utilities/boosthead.hpp>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>


namespace util {
  namespace num {
    namespace lapacke {



      /*
       * SYEV - solve the eigensystem for a real symmetric matrix
       */
      template <typename T>
      size_t syev(const ublas::matrix<T> & mat, ublas::vector<T> & eval, ublas::matrix<T> & evec);

      template <>
      inline
      size_t syev(const ublas::matrix<double> & mat, ublas::vector<double> & eval, ublas::matrix<double> & evec);

      template <>
      inline
      size_t syev(const ublas::matrix<float> & mat, ublas::vector<float> & eval, ublas::matrix<float> & evec);



//----------------------------------------------------------------------------------
//----------------implementations---------------------------------------------------
//----------------------------------------------------------------------------------



      template <>
      //inline
      size_t syev(const ublas::matrix<double> & mat, ublas::vector<double> & eval, ublas::matrix<double> & evec){




        //evec.resize(mat.size1(), mat.size2());
        evec = mat;
        eval.resize(evec.size1());

        size_t res = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', evec.size1(), &evec.data()[0], evec.size1(), &eval.data()[0]);
        return res;

      }

      template <>
      //inline
      size_t syev(const ublas::matrix<float> & mat, ublas::vector<float> & eval, ublas::matrix<float> & evec){
        evec = mat;
        eval.resize(evec.size1());
        size_t res = LAPACKE_ssyev( LAPACK_ROW_MAJOR, 'V', 'U', evec.size1(), &evec.data()[0], evec.size1(), &eval.data()[0]);
        return res;

      }





      }
  }
}






#endif // LAPACK_HPP
