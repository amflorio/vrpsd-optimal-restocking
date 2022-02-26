#ifndef ROUNDEDCAPACITYCUT_H
#define ROUNDEDCAPACITYCUT_H

#include <vector>
#include "Arc.h"

class RoundedCapacityCut {
    public:
        std::vector<Arc> arcs;
        int rhs;
        RoundedCapacityCut(std::vector<Arc> cutarcs, int cutrhs) :
                arcs{move(cutarcs)}, rhs{cutrhs} {}
        void print() const;
};

#endif

