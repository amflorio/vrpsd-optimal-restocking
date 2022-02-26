#include "RoundedCapacityCut.h"

using namespace std;

void RoundedCapacityCut::print() const {
    cout<<"arcs: ";
    for (const auto& a : arcs)
        cout<<a<<", ";
    cout<<"\trhs: "<<rhs<<endl;
}

