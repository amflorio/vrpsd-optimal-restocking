#include <cmath>
#include <iostream>
#include "Binomial.h"
#include "NegBinomial.h"

using namespace std;

// See NegBin.pdf file for the parametrization used
double NegBinomial::pmf(int x, int mean, double p) {
    int r;
    if (p==0.5)
        r=mean;
    else if (p==2.0/3.0)
        r=2*mean;
    else {
        cerr<<"NegBinomial::pmf(): invalid p"<<endl;
        exit(-1);
    }
    return Binomial::nchoosek(x+r-1, x)*pow(p, r)*pow(1-p, x);
}

