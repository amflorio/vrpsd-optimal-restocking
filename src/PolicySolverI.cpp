#include <cmath>
#include <iostream>
#include "Instance.h"
#include "PolicySolverI.h"
#include "Route.h"

using namespace std;

PolicySolverI::PolicySolverI(const Instance& i, const Route& r) : inst{i},
        apriori{r} {
}

double PolicySolverI::solve() const {
    vector<int> nodes=apriori.nodes();
    if (nodes.size()<3 || nodes.at(0)!=0 || nodes.back()!=0) {
        cerr<<"attempted to compute SDP on invalid route"<<endl;
        cerr<<"route: "<<apriori<<endl;
        exit(-1);
    }
    vector<double> vq;
    int pos=nodes.size()-2;
    int Q=inst.vehicleQ();
    // base case
    for (int q=0; q<=Q; ++q)
        vq.push_back(inst.dist(nodes.at(pos), 0));
    pos--;
    // general case
    while (pos>0) {
        vector<double> newvq;
        // indirect (replenish first at the depot)
        double depot=inst.dist(nodes.at(pos), 0)+inst.dist(0, nodes.at(pos+1));
        for (int k=0; k<=inst.maxDemand(nodes.at(pos+1)); ++k) {
            depot+=inst.prob(nodes.at(pos+1), k)*((inst.dist(nodes.at(pos+1), 0)
                    +inst.dist(0, nodes.at(pos+1)))*gamma(k, Q)
                    +vq.at(Q*(gamma(k, Q)+1)-k));
        }
        for (int q=0; q<=Q; ++q) {
            // direct
            double direct=inst.dist(nodes.at(pos), nodes.at(pos+1));
            for (int k=0; k<=inst.maxDemand(nodes.at(pos+1)); ++k) {
                direct+=inst.prob(nodes.at(pos+1), k)*
                        ((inst.dist(nodes.at(pos+1), 0)
                        +inst.dist(0, nodes.at(pos+1)))*gamma(k, q)
                        +vq.at(Q*gamma(k, q)+q-k));
            }
            newvq.push_back(std::min(direct, depot));
        }
        vq=newvq;
        pos--;
    }
    // final case
    double ecost=inst.dist(0, nodes.at(1));
    for (int k=0; k<=inst.maxDemand(nodes.at(1)); ++k) {
        ecost+=inst.prob(nodes.at(1), k)*((inst.dist(nodes.at(1), 0)
                +inst.dist(0, nodes.at(1)))*gamma(k, Q)
                +vq.at(Q*(gamma(k, Q)+1)-k));
    }
    return ecost;
}

