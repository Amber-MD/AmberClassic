#ifndef timer_struc_HPP
#define timer_struc_HPP

#include <iostream>
using namespace std;

namespace rism3d_c{

    struct timer_cpp{
        string* name;
        bool running;
        bool startedParent;
        bool displaced;
        int sublevels;
        double total;
        double timestamp;
        timer_cpp *parent;
        timer_cpp **child;
    };
    typedef struct timer_cpp timer;
}

#endif