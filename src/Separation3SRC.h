#ifndef SEPARATION3SRC_H
#define SEPARATION3SRC_H

#include <vector>
#include "Solution.h"
#include "SubsetRowCut.h"

class Separation3SRC {
    private:
        static std::shared_ptr<Data> data;
        Solution solution;
    public:
        Separation3SRC(Solution sol) : solution{std::move(sol)} {}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        std::vector<SubsetRowCut> solve() const;
};

#endif

