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


#ifndef STATISTICS_HPP
#define STATISTICS_HPP


namespace utils {
  namespace stats{


  /*
   * Class for computing statistics in one pass
   *
   */
    template <typename T>
    class SPStats
    {
    public:
      // ctor
      SPStats();

      // start over
      void clear();

      // add each data point
      void push(T x);

      // add data in full
      template<typename TD> void pushall(TD & data);

      // get the size of the accumulated data
      size_t size() const;

      // mean
      T mean() const;

      // variance
      T variance() const;

      // stddev
      T standarddeviation() const;

      // stderr
      T standarderror() const;

      // skewness
      T skewness() const;

      // kurtoisis
      T kurtosis() const;

    private:
      size_t n;
      T M1, M2, M3, M4;
    };


template <typename T>
  SPStats<T>::SPStats()
  {
      clear();
  }

template <typename T>
  void SPStats<T>::clear()
  {
      n = 0;
      M1 = M2 = M3 = M4 = 0.0;
  }

template <typename T>
  void SPStats<T>::push(T x)
  {
      T delta, delta_n, delta_n2, term1;

      size_t n1 = n;
      n++;
      delta = x - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * n1;
      M1 += delta_n;
      M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
  }

template <typename T>
  size_t SPStats<T>::size() const
  {
      return n;
  }

template <typename T>
  T SPStats<T>::mean() const
  {
      return M1;
  }

template <typename T>
  T SPStats<T>::variance() const
  {
      return M2/(n-1.0);
  }

template <typename T>
  T SPStats<T>::standarddeviation() const
  {
      return std::sqrt( this->variance() );
  }


  template <typename T>
    T SPStats<T>::standarderror() const
    {
        return (this->standarddeviation())/std::sqrt(T(n)) ;
    }



template <typename T>
  T SPStats<T>::skewness() const
  {
      return sqrt(T(n)) * M3/ pow(M2, 1.5);
  }

template <typename T>
  T SPStats<T>::kurtosis() const
  {
      return T(n)*M4 / (M2*M2) - 3.0;
  }



template<typename T>
template<typename TD>
  void SPStats<T>::pushall(TD & data){
    for (size_t i = 0 ; i != data.size(); ++i){
        this->push(data[i]);
      }
  }


  }
}




#endif // STATISTICS_HPP
