#ifndef ARC_H
#define ARC_H

#include <iostream>

class Arc {
    public:
        int i;
        int j;
        double val;
        Arc(int a_i, int a_j, double v) : i{a_i}, j{a_j}, val{v} {}
};

std::ostream& operator<<(std::ostream& os, const Arc& a);

#endif

