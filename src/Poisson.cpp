#include <cmath>
#include <iostream>
#include "Poisson.h"

using namespace std;

double Poisson::pmf(const double k, const double lambda) {
    double ret=exp(k*log(lambda)-lgamma(k+1.0)-lambda);
    if (std::isinf(ret) || std::isnan(ret) || ret<0 || ret>1) {
        cerr<<"pmf("<<k<<", "<<lambda<<"): invalid prob.: "<<ret<<endl;
        exit(-1);
    }
    return ret;
}

