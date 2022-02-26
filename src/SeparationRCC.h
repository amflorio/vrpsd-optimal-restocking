#ifndef SEPARATIONRCC_H
#define SEPARATIONRCC_H

#include <memory>
#include "Data.h"
#include "LinearSolver.h"
#include "Solution.h"

class SeparationRCC {
    private:
        static std::shared_ptr<Data> data;
        Solution sol;
        void checkCut(const RoundedCapacityCut& cut,
                const std::vector<int>& nodes) const;
    public:
        SeparationRCC(const Solution& s) : sol{s} {}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        std::vector<RoundedCapacityCut> solve() const;
};

#endif

