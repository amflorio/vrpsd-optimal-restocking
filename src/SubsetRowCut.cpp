#include <cassert>
#include "SubsetRowCut.h"

using namespace std;

shared_ptr<Data> SubsetRowCut::data;

SubsetRowCut::SubsetRowCut(vector<int> set, double rhs, double v)
        : set_{move(set)}, contains_(data->numNodes(), 0), rhs_{rhs}, vio{v} {
    for (auto i : set_)
        contains_[i]=1;
}

ostream& operator<<(ostream& os, const SubsetRowCut& c) {
    assert(c.set_.size()>0);
    os<<"[{"<<c.set_[0];
    for (int i=1; i<c.set_.size(); ++i)
        os<<","<<c.set_[i];
    os<<"} rhs="<<c.rhs_<<" vio="<<c.vio<<"]";
    return os;
}

