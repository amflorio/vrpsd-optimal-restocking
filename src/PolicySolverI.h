#ifndef POLICYSOLVERI_H
#define POLICYSOLVERI_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "Instance.h"
#include "Route.h"

class PolicySolverI {
    private:
        Instance inst;
        Route apriori;
        int gamma(int k, int q) const {return std::max(0,
                1+(int)std::floor(((k-q-1)*1.0)/inst.vehicleQ()));}
    public:
        PolicySolverI(const Instance& i, const Route& r);
        double solve() const;
};

#endif

