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

#ifndef READDX_HPP_
#define READDX_HPP_

#include "boosthead.hpp"
#include <fstream>
#include <iostream>


namespace util {
    
    namespace files {
        
        
        template <typename T, typename Tstream>
        void readxyzvgen(Tstream & stream, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                         boost::multi_array<T, 3> &data3d) {
            
            
            orig.resize(3);
            rot.resize(3, 3);
            
            orig *= 0.0;
            rot *= 0.0;
            
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            typedef tokenizer::iterator tokenizerit;
            boost::char_separator<char> sep(" ");
            
            
            std::string line;
            
            
            size_t linecounter = 1;
            //   util::convto<size_t> tst;
            //   util::convto<float> tft;
            
            size_t /*n = 0,*/ nmax = 0;
            
            while (getline(stream, line)) {
                
                tokenizer tok(line, sep);
                
                if (linecounter > 3) {
                    tokenizerit beg = tok.begin();
                    try {// just in case straneous things are written in the file.
                        size_t n1 = boost::lexical_cast<size_t>(*beg);
                        ++beg;
                        size_t n2 = boost::lexical_cast<size_t>(*beg);
                        ++beg;
                        size_t n3 = boost::lexical_cast<size_t>(*beg);
                        ++beg;
                        data3d[n1][n2][n3] += boost::lexical_cast<T>(*beg);
                    } catch (...) {};
                } else if (linecounter == 1) {
                    // here we read the data dimensions
                    tokenizerit beg = tok.begin();
                    ++beg;
                    size_t n1 = boost::lexical_cast<size_t>(*beg);
                    ++beg;
                    size_t n2 = boost::lexical_cast<size_t>(*beg);
                    ++beg;
                    size_t n3 = boost::lexical_cast<size_t>(*beg);
                    nmax = n1 * n2 * n3;
                    
                    if (data3d.num_elements()==0) {
                        data3d.resize(boost::extents[n1][n2][n3]);
                        for (size_t i = 0; i != data3d.num_elements(); ++i) data3d.data()[i] = 0.0;
                    } else if (data3d.num_elements() != nmax) {
                        std::exit(0);
                    } else {
                        std::cout << "# accumulating data" << std::endl;
                    }
                } else if (linecounter == 2) {
                    // here we read the origin
                    tokenizerit beg = tok.begin();
                    ++beg;
                    orig[0] = boost::lexical_cast<T>(*beg);
                    ++beg;
                    orig[1] = boost::lexical_cast<T>(*beg);
                    ++beg;
                    orig[2] = boost::lexical_cast<T>(*beg);
                } else if (linecounter == 3) {
                    // here we read the rotation matrix
                    tokenizerit beg = tok.begin();
                    ++beg;
                    rot(0, 0) = boost::lexical_cast<T>(*beg);
                    ++beg;
                    rot(1, 1) = boost::lexical_cast<T>(*beg);
                    ++beg;
                    rot(2, 2) = boost::lexical_cast<T>(*beg);
                }
                
                
                linecounter++;
            }
            
            
        }
        
        
        
        template<typename T, typename Tstream>
        
        void
        readdxgen(Tstream &stream, ublas::vector <T> &orig, ublas::matrix <T> &rot, boost::multi_array<T, 3> &data3d) {
            
            
            orig.resize(3);
            rot.resize(3, 3);
            
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            typedef tokenizer::iterator tokenizerit;
            boost::char_separator<char> sep(" ");
            
            typedef std::vector<std::string> Tvs;
            Tvs stokens;
            
            std::string line;
            size_t linecounter = 1;
            util::convto<size_t> tst;
            util::convto<float> tft;
            
            
            size_t n = 0, nmax = 0;
            
            
            //        bool isalloc = false;
            
            
            while (getline(stream, line)) {
                if (line[0] != '#') {
                    if (linecounter < 8) {
                        // here we read the header
                        stokens = boost::copy_range<Tvs>(
                                                         line | boost::adaptors::tokenized(std::string("[\\t\\s,]+"), -1));
                        if (linecounter == 1) {
                            size_t n1 = tst(stokens[5]);
                            size_t n2 = tst(stokens[6]);
                            size_t n3 = tst(stokens[7]);
                            nmax = n1 * n2 * n3;
                            
                            if (data3d.num_elements() == 0) {
                                data3d.resize(boost::extents[n1][n2][n3]);
                                for (size_t i = 0; i != data3d.num_elements(); ++i) data3d.data()[i] = 0.0;
                            } else if (data3d.num_elements() != nmax) {
                                std::exit(0);
                            } else {
                                std::cout << "# accumulating data" << std::endl;
                            }
                       } else if (linecounter == 2) {
                            orig[0] = tft(stokens[1]);
                            orig[1] = tft(stokens[2]);
                            orig[2] = tft(stokens[3]);
                       /* } else if (linecounter == 3) {
                            rot(0, 0) = tft(stokens[1]);
                            rot(0, 1) = tft(stokens[2]);
                            rot(0, 2) = tft(stokens[3]);
                        } else if (linecounter == 4) {
                            rot(1, 0) = tft(stokens[1]);
                            rot(1, 1) = tft(stokens[2]);
                            rot(1, 2) = tft(stokens[3]);
                        } else if (linecounter == 5) {
                            rot(2, 0) = tft(stokens[1]);
                            rot(2, 1) = tft(stokens[2]);
                            rot(2, 2) = tft(stokens[3]);
                        }*/
                        
                        
                    } else if (linecounter == 3) {
                        rot(0, 0) = tft(stokens[1]);
                        rot(1, 0) = tft(stokens[2]);
                        rot(2, 0) = tft(stokens[3]);
                    } else if (linecounter == 4) {
                        rot(0, 1) = tft(stokens[1]);
                        rot(1, 1) = tft(stokens[2]);
                        rot(2, 1) = tft(stokens[3]);
                    } else if (linecounter == 5) {
                        rot(0, 2) = tft(stokens[1]);
                        rot(1, 2) = tft(stokens[2]);
                        rot(2, 2) = tft(stokens[3]);
                    }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    } else {
                        if (n < nmax) {
                            
                            tokenizer tok(line, sep);
                            for (tokenizerit beg = tok.begin(); beg != tok.end(); ++beg) {
                                try {
                                    data3d.data()[n] += boost::lexical_cast<T>(*beg);
                                } catch (const std::exception& e) {
                                    
                                    //std::cout << line << " "<< *beg << std::endl;
                                }
                                
                                if (data3d.data()[n] != data3d.data()[n]) {
                                    
                                    std::cout << "# non numeric values!" << std::endl;
                                    std::cout << linecounter << " " << line << std::endl;
                                    std::exit(0);
                                };
                                
                                ++n;
                            }
                        }
                    }
                    linecounter++;
                }
            }
            
            
        };
        
        
        template<typename T>
        bool readdx(const std::string &filename, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                    boost::multi_array<T, 3> &data3d) {
            
            
            std::ifstream file(filename.c_str());
            
            if (file) {
                //std::cout << "# Reading from " << filename << " file ." << std::endl;
                
                readdxgen(file, orig, rot, data3d);
                
                file.close();
                
                return true;
            } else {
                
                std::cout << ">> Cannot open:  " << filename << std::endl;
                return false;
            }
            
            
        };
        
        
        template<typename T>
        bool readdxbz2(const std::string &filename, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                       boost::multi_array<T, 3> &data3d) {
            
            
            std::ifstream file(filename.c_str(), std::ios_base::in | std::ios_base::binary);
            
            
            if (file) {
                std::cout << "# Reading from " << filename << " file ." << std::endl;
                boost::iostreams::filtering_istream in;
                in.push(boost::iostreams::bzip2_decompressor());
                in.push(file);
                
                readdxgen(in, orig, rot, data3d);
                
                file.close();
                
                return true;
            } else {
                
                std::cout << ">> Cannot open:  " << filename << std::endl;
                return false;
            }
            
            
        };
        
        
        template<typename T>
        bool readdxgz(const std::string &filename, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                      boost::multi_array<T, 3> &data3d) {
            
            
            std::ifstream file(filename.c_str(), std::ios_base::in | std::ios_base::binary);
            
            
            if (file) {
                std::cout << "# Reading from " << filename << " file ." << std::endl;
                boost::iostreams::filtering_istream in;
                in.push(boost::iostreams::gzip_decompressor());
                in.push(file);
                
                
                readdxgen(in, orig, rot, data3d);
                
                file.close();
                
                return true;
            } else {
                
                std::cout << ">> Cannot open:  " << filename << std::endl;
                return false;
            }
            
            
        };
        
        
        template<typename T>
        bool readxyzv(const std::string &filename, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                      boost::multi_array<T, 3> &data3d) {
            
            std::ifstream file(filename.c_str());
            
            if (file) {
                std::cout << "# Reading from " << filename << " file." << std::endl;
                
                readxyzvgen(file,orig,rot,data3d);
                
                file.close();
                return true;
            } else {
                
                std::cout << ">> Cannot open:  " << filename << std::endl;
                return false;
            }
            
            
        }
        
        
        template<typename T>
        bool readxyzvgz(const std::string &filename, ublas::vector <T> &orig, ublas::matrix <T> &rot,
                        boost::multi_array<T, 3> &data3d) {
            
            std::ifstream file(filename.c_str(), std::ios_base::in | std::ios_base::binary);
            
            if (file) {
                std::cout << "# Reading from " << filename << " file." << std::endl;
                
                boost::iostreams::filtering_istream in;
                in.push(boost::iostreams::gzip_decompressor());
                in.push(file);
                
                readxyzvgen(in,orig,rot,data3d);
                
                file.close();
                return true;
            } else {
                
                std::cout << ">> Cannot open:  " << filename << std::endl;
                return false;
            }
            
            
        }
        
        
    }
}

#endif /* READDX_HPP_ */

