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



#ifndef ROTATIONS_HPP
#define ROTATIONS_HPP


#include <iostream>
#include <boost/math/quaternion.hpp>
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>
namespace ublas = boost::numeric::ublas;

namespace rotations{


    void rotmatrixaxis(double theta, ublas::vector<double> & axis, ublas::matrix<double>  & rotm){
        
        rotm.resize(3,3);
        ublas::matrix<double>  ux(3,3,0), utens(3,3,0);
        
        ux(0,1) = -axis(2);
        ux(0,2) =  axis(1);
        ux(1,0) =  axis(2);
        ux(1,2) = -axis(0);
        ux(2,0) = -axis(1);
        ux(2,1) =  axis(0);
        
        utens = ublas::outer_prod(axis,axis);
        ublas::identity_matrix<double> id(3);
        rotm = id * std::cos(theta) + ux * std::sin(theta) + (1-std::cos(theta)) * utens;
        
    }
    
    
    
  template <typename T>
  ublas::vector<T>
  rotatevector(const ublas::vector<T> & v, const boost::math::quaternion<T> & q ){

      boost::math::quaternion<T> qv (0,v(0),v(1),v(2));
      ublas::vector<T> vr(3);
      qv  =  q * qv * boost::math::conj(q) ;
      vr(0) = qv.R_component_2();
      vr(1) = qv.R_component_3();
      vr(2) = qv.R_component_4();
      return vr;

  }


  template <typename T>
  void slerp(const boost::math::quaternion<T> & q1,
             const boost::math::quaternion<T> & q2,
                   T alpha,
                   boost::math::quaternion<T> & q
             ){

   T beta; // complementary interp parameter
   T theta; // angle between quaternion
   T sin_t,cos_t; // sine and cosine
   bool bflip = false;

   cos_t =
       q1.R_component_1() * q2.R_component_1() +
       q1.R_component_2() * q2.R_component_2() +
       q1.R_component_3() * q2.R_component_3() +
       q1.R_component_4() * q2.R_component_4()
       ;


   if (cos_t < 0){
       cos_t *= -1.0;
       bflip=true;
     }

   if (1.0 - cos_t < 1e-7 ){

       beta =  1.0 - alpha;

     } else {

       theta = std::acos(cos_t);
       sin_t = std::sin(theta);
       beta  = std::sin(theta - alpha*theta) / sin_t;
       alpha = std::sin(alpha*theta) / sin_t;

     }

   if (bflip) alpha *= -1.0;

   q = beta*q1+alpha*q2;

   //q /= boost::math::abs(q);
  }

template <typename T>
void CatRom(
    const boost::math::quaternion<T> & q00,
    const boost::math::quaternion<T> & q01,
    const boost::math::quaternion<T> & q02,
    const boost::math::quaternion<T> & q03,
          boost::math::quaternion<T> & q,
          T t
    ){

 boost::math::quaternion<T>  q10,q11,q12,q20,q21;
 slerp(q00,q01,t+1,q10);
 slerp(q01,q02,t,q11);
 slerp(q02,q03,t-1,q12);

 slerp(q10,q11,(t+1)/2.0,q20);
 slerp(q11,q12,t/2.0,q21);
 slerp(q20,q21,t,q);
}

template <typename T>
void UBS(
    const boost::math::quaternion<T> & q00,
    const boost::math::quaternion<T> & q01,
    const boost::math::quaternion<T> & q02,
    const boost::math::quaternion<T> & q03,
          boost::math::quaternion<T> & q,
          T t
    ){

 boost::math::quaternion<T>  q10,q11,q12,q20,q21;
 slerp(q00,q01,(t+2)/3,q10);
 slerp(q01,q02,(t+1)/3,q11);
 slerp(q02,q03,t/3,q12);

 slerp(q10,q11,(t+1)/2.0,q20);
 slerp(q11,q12,t/2.0,q21);
 slerp(q20,q21,t,q);
}








template <typename T>
boost::math::quaternion<T>
log(const boost::math::quaternion<T> & q ){
  /*
   *
   * log of quaternion as implemented in opengl
   *
   *
  tvec3<T, P> u(q.x, q.y, q.z);
  T Vec3Len = length(u);

  if (Vec3Len < epsilon<T>())
  {
          if(q.w > static_cast<T>(0))
                  return tquat<T, P>(log(q.w), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
          else if(q.w < static_cast<T>(0))
                  return tquat<T, P>(log(-q.w), pi<T>(), static_cast<T>(0), static_cast<T>(0));
          else
                  return tquat<T, P>(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
  }
  else
  {
          T QuatLen = sqrt(Vec3Len * Vec3Len + q.w * q.w);
          T t = atan(Vec3Len, T(q.w)) / Vec3Len;
          return tquat<T, P>(log(QuatLen), t * q.x, t * q.y, t * q.z);
  }

  */


  T complength = std::sqrt(  q.R_component_2()* q.R_component_2() +
                             q.R_component_3()* q.R_component_3() +
                             q.R_component_4()* q.R_component_4()
                             );

   if (complength < 1e-7){

       if(q.R_component_1() > 0.0){
             return boost::math::quaternion<T>(std::log(q.R_component_1()),0.0,0.0,0.0);
         }
         else if(q.R_component_1() > 0.0){
             return boost::math::quaternion<T>(std::log(-q.R_component_1()),M_PI,0.0,0.0);
         }
         else{
             return boost::math::quaternion<T>(std::numeric_limits<T>::infinity(),
                                             std::numeric_limits<T>::infinity(),
                                             std::numeric_limits<T>::infinity(),
                                             std::numeric_limits<T>::infinity()
                                             );
         }
     }
   else {

       T t  = std::atan2(complength,q.R_component_1())/complength;
       //std::cout << 	boost::math::abs(q) << std::endl;
       return boost::math::quaternion<T>( std::log(boost::math::abs(q)),
                                          t*q.R_component_2(),
                                          t*q.R_component_3(),
                                          t*q.R_component_4()
                                          );

     }


}

template <typename T>
void CRspline(
    const boost::math::quaternion<T> & q0,
    const boost::math::quaternion<T> & q1,
    const boost::math::quaternion<T> & q2,
    const boost::math::quaternion<T> & q3,
          boost::math::quaternion<T> & q,
          T t
    ){

  boost::math::quaternion<T> w1,w2,w3;

  w1 = rotations::log(q0*boost::math::conj(q2));
  w3 = rotations::log(q1*boost::math::conj(q3));
  //w2 = rotations::log(q1*boost::math::conj(q2));
  w1 /= 3.0;
  w3 /= 3.0;

  w2=rotations::log (
        boost::math::conj(boost::math::exp(w1)) *
        boost::math::conj(q1) *
        q2 *
        boost::math::conj(boost::math::exp(w3))
        );

  q = q1 * boost::math::exp(w1*(1-(1-t)*(1-t)*(1-t))) *
      boost::math::exp(w2*(3*t*t-2*t*t*t))*
      boost::math::exp(w3*t*t*t);

}


/*
  typedef boost::random::mt19937 Tmt;
  Tmt rng(std::time(NULL)), rngo(std::time(NULL));

  template <typename T>
  ublas::vector<T> randtrans(T & incr){
          boost::random::normal_distribution<> xran(0,incr);
          T r = xran(rng);
          typedef ublas::vector<T> TVec;
          typedef boost::uniform_on_sphere<T, TVec > Tusph;
          boost::variate_generator<Tmt,Tusph > vgen( rngo,Tusph(3));
          TVec s = vgen();
          s *= r;
          return s;

  }
*/

/*
  template<typename T>
  boost::math::quaternion<T> randquat(T & rdisp){
          // choose sphere radius
          boost::random::normal_distribution<> xran(0,rdisp);
          T r = xran(rng);

          // choose vector sphere
          typedef ublas::vector<T> TVec;
          typedef boost::uniform_on_sphere<T, TVec > Tusph;
          boost::variate_generator<Tmt,Tusph > vgen( rngo,Tusph(3));

          TVec s = vgen();
          s *= r;

          // return exp

          boost::math::quaternion<T> q(0.0,s[0],s[1],s[2]);
          return boost::math::exp(q);


  }
*/
  template<typename T>
  ublas::matrix<T> qtorot(const boost::math::quaternion<T> & q){

          boost::math::quaternion<T> qt(q);

          qt /= boost::math::abs(qt);


          ublas::matrix<T> m(3,3,0.0);
          T qw = qt.R_component_1();
          T qx = qt.R_component_2();
          T qy = qt.R_component_3();
          T qz = qt.R_component_4();

          m(0,0) = 1.0f - 2.0f*qy*qy - 2.0f*qz*qz ;
          m(0,1) = 2.0f*qx*qy - 2.0f*qz*qw;
          m(0,2) = 2.0f*qx*qz + 2.0f*qy*qw;
          //m(0,3) = 0.0f;
          m(1,0) = 2.0f*qx*qy + 2.0f*qz*qw;
          m(1,1) = 1.0f - 2.0f*qx*qx - 2.0f*qz*qz ;
          m(1,2) = 2.0f*qy*qz - 2.0f*qx*qw;
          //m(1,3) = 0.0f;
          m(2,0) = 2.0f*qx*qz - 2.0f*qy*qw;
          m(2,1) = 2.0f*qy*qz + 2.0f*qx*qw;
          m(2,2) = 1.0f - 2.0f*qx*qx - 2.0f*qy*qy ;
          //m(2,3) = 0.0f;
         // m(3,0) = 0.0f;
         // m(3,1) = 0.0f;
         // m(3,2) = 0.0f;
         // m(3,3) = 1.0f;


          return m;
  }



  template <typename T>
  ublas::matrix<T> ttorott(ublas::vector<T> & v){
          ublas::matrix<T> m(4,4,0.0);
          m(3,0) = v(0);
          m(3,1) = v(1);
          m(3,2) = v(2);
          m(3,3) = 1.0;
          return m;

  }

  /*
   * Convert a rotation matrix to a
   * quaternion.
   */

  template <typename T>
  boost::math::quaternion<T>
  convrotmatquat(const ublas::matrix<T> & m)
  {

    T tr = m(0,0) + m(1,1) + m(2,2);
    T qw,qx,qy,qz,S;

    if (tr > 0) {
      S = std::sqrt(tr+1.0) * 2; // S=4*qw
      qw = 0.25 * S;
      qx = (m(2,1) - m(1,2)) / S;
      qy = (m(0,2) - m(2,0)) / S;
      qz = (m(1,0) - m(0,1)) / S;
    } else
    if ((m(0,0) > m(1,1))&(m(0,0) > m(2,2)))
    {
      S = std::sqrt(1.0 + m(0,0) - m(1,1) - m(2,2)) * 2; // S=4*qx
      qw = (m(2,1) - m(1,2)) / S;
      qx = 0.25 * S;
      qy = (m(0,1) + m(1,0)) / S;
      qz = (m(0,2) + m(2,0)) / S;
    } else
    if (m(1,1) > m(2,2))
    {
      S = std::sqrt(1.0 + m(1,1) - m(0,0) - m(2,2)) * 2; // S=4*qy
      qw = (m(0,2) - m(2,0)) / S;
      qx = (m(0,1) + m(1,0)) / S;
      qy = 0.25 * S;
      qz = (m(1,2) + m(2,1)) / S;
    } else
    {
      S = std::sqrt(1.0 + m(2,2) - m(0,0) - m(1,1)) * 2; // S=4*qz
      qw = (m(1,0) - m(0,1)) / S;
      qx = (m(0,2) + m(2,0)) / S;
      qy = (m(1,2) + m(2,1)) / S;
      qz = 0.25 * S;
    }

      return boost::math::quaternion<T>(qw,qx,qy,qz);

  }


  template <typename T>
  void convrotmatquat(
      const std::vector< ublas::matrix<T> > & mat,
      ublas::vector< boost::math::quaternion<T> > & q
  )
  {
    q.resize(mat.size());

     for (size_t i = 0 ; i != mat.size();++i){
         q(i) = convrotmatquat(mat[i]);
     }

  }


  template <typename T>
  void convquatrotmat(
      const ublas::vector< boost::math::quaternion<T> > & q,
      ublas::vector< ublas::matrix<T> > & mat

  )
  {
    mat.resize(q.size());

     for (size_t i = 0 ; i != mat.size();++i){
         mat(i) = qtorot(q(i));
     }

  }





}


#endif // ROTATIONS_HPP
