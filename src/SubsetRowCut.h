#ifndef SUBSETROWCUT_H
#define SUBSETROWCUT_H

#include <iostream>
#include <vector>
#include "Data.h"

class SubsetRowCut {
    friend std::ostream& operator<<(std::ostream& os, const SubsetRowCut& c);
    private:
        static std::shared_ptr<Data> data;
        std::vector<int> set_;
        std::vector<int> contains_;
        double rhs_;
        double vio;
    public:
        SubsetRowCut(std::vector<int> set, double rhs, double v);
        bool contains(int i) const {return contains_[i];}
        std::vector<int> nodes() const {return set_;}
        double rhs() const {return rhs_;}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        double violation() const {return vio;}
};

#endif

std::ostream& operator<<(std::ostream& os, const SubsetRowCut& c);

