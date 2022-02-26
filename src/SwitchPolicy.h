#ifndef SWITCHPOLICY_H
#define SWITCHPOLICY_H

#include <vector>
#include "Route.h"

class SwitchPolicy {
    private:
        static std::shared_ptr<Data> data;
        Route apriori;
        double eLength=-1;
        std::vector<std::vector<std::vector<double>>> nu;
        void computePolicy();
        int gamma(int k, int q) const {return std::max(0,
                1+(int)std::floor(((k-q-1)*1.0)/data->vehicleQ()));}
        double piprime(int i, int j, int f, int q) const;
        double pi2prime(int i, int j, int f) const;
        double pistar(int i, int j, int f, int q) const
                {return std::min(piprime(i,j,f,q),pi2prime(i,j,f));}
    public:
        SwitchPolicy(Route r);
        double solve() {return eLength;}
        static void setData(std::shared_ptr<Data> d) {data=d;}
};

#endif

