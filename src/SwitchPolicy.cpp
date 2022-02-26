#include <cassert>
#include <cmath>
#include <iostream>
#include "Data.h"
#include "Route.h"
#include "SwitchPolicy.h"

using namespace std;

shared_ptr<Data> SwitchPolicy::data;

SwitchPolicy::SwitchPolicy(Route r) : apriori{move(r)} {
    assert(data!=nullptr);
    //assert(apriori.numVisits()>1);
    computePolicy();
    const int Q=data->vehicleQ();
    if (apriori.numVisits()>1)
        eLength=min(piprime(0,apriori[1],1,Q), piprime(0,apriori[2],1,Q));
    else
        eLength=piprime(0,apriori[1],1,Q);
}

void SwitchPolicy::computePolicy() {
    const int s_fminus1=0;
    const int s_f=1;
    const int s_fplus1=2;
    const int F=apriori.numVisits();
    // data structure initialization - nu(n,F,q)
    nu=vector<vector<vector<double>>>(3, vector<vector<double>>(F+1,
            vector<double>(data->vehicleQ()+1, -1.0)));
    // data structure initialization - nu_decisions(n,F,q) - SKIPPING
    // dynamic programming base case
    for (int q=0; q<=data->vehicleQ(); ++q) {
        nu[s_fminus1][F][q]=data->dist(apriori[F-1], 0);
        nu[s_f][F][q]=data->dist(apriori[F], 0);
    }
    // dynamic programming general case
    for (int f=F-1; f>=1; --f) {
        for (int q=0; q<=data->vehicleQ(); ++q) {
            // if n=s_{f+1}
            nu[s_fplus1][f][q]=pistar(apriori[f+1],apriori[f],f+1,q);
            // if f=F-1
            if (f==F-1) {
                nu[s_fminus1][f][q]=pistar(apriori[f-1],apriori[F],F,q);
                nu[s_f][f][q]=pistar(apriori[f],apriori[F],F,q);
            } else {
                nu[s_fminus1][f][q]=min(pistar(apriori[f-1],apriori[f+1],f+1,q),
                        pistar(apriori[f-1],apriori[f+2],f+1,q));
                nu[s_f][f][q]=min(pistar(apriori[f],apriori[f+1],f+1,q),
                        pistar(apriori[f],apriori[f+2],f+1,q));
            }
        }
    }
}

double SwitchPolicy::piprime(int i, int j, int f, int q) const {
    double cij=data->dist(i, j);
    double c0j=data->dist(0, j);
    double cj0=data->dist(j, 0);
    int Q=data->vehicleQ();
    assert(j==apriori[f-1] || j==apriori[f] || j==apriori[f+1]);
    int idxj=j==apriori[f-1]?0:j==apriori[f]?1:2;
    double cost=cij;
    for (int k=0; k<=data->maxDemand(j); ++k) {
        int gkq=gamma(k, q);
        assert(nu.at(idxj).at(f).at(q+Q*gkq-k)>=0);
        cost+=(gkq*(cj0+c0j)+nu[idxj][f][q+Q*gkq-k])*data->prob(j, k);
    }
    return cost;
}

double SwitchPolicy::pi2prime(int i, int j, int f) const {
    double ci0=data->dist(i, 0);
    double c0j=data->dist(0, j);
    double cj0=data->dist(j, 0);
    int Q=data->vehicleQ();
    assert(j==apriori[f-1] || j==apriori[f] || j==apriori[f+1]);
    int idxj=j==apriori[f-1]?0:j==apriori[f]?1:2;
    double cost=ci0+c0j;
    for (int k=0; k<=data->maxDemand(j); ++k) {
        int gkQ=gamma(k, Q);
        assert(nu.at(idxj).at(f).at(Q+Q*gkQ-k)>=0);
        cost+=(gkQ*(cj0+c0j)+nu[idxj][f][Q+Q*gkQ-k])*data->prob(j, k);
    }
    return cost;
}

