#include <cmath>
#include "Binomial.h"

using namespace std;

double Binomial::nchoosek(int n, int k) {
    double ret=1.0;
    for (int i=1; i<=k && i<=n-k; ++i)
        ret*=(n+1-i)/(i*1.0);
    return ret;
}

double Binomial::pmf(int k, int mean) {
    // mean=np
    double p=0.5;
    // 2 hardcoded. Other places hardcoded as well (e.g. VRPSDData, Instance).
    int n=2*mean;
    return nchoosek(n, k)*pow(p, k)*pow(1-p, n-k);
}

