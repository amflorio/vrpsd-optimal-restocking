#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include "Config.h"
#include "Data.h"
#include "Label.h"
#include "PolicySolverI.h"

using namespace std;

Data::Data(const Instance& i) : instance{i} {
    genDemands();
    genMinDists();
    genMinDistsPlus();
    // approximate size of a Label (in bytes)
    sizeRLabel=sizeof(Label);
    sizeULabel=sizeRLabel;
    sizeULabel+=numNodes()*sizeof(int);             // _visitsNode
    sizeULabel+=10*sizeof(int);                     // estimate: routeNodes
    sizeULabel+=vehicleQ()*sizeof(double);          // vq
    cout<<"sizeReservedLabel="<<sizeRLabel<<endl;
    cout<<"sizeUsedLabel="<<sizeULabel<<endl;
}

int Data::countBits(int s) const {
    int count=0;
    while (s>0) {
        if ((s&1)==1)
            count++;
        s>>=1;
    }
    return count;
}

/* Organize the wrappers around the demand information from Instance */
void Data::genDemands() {
    const double MY_ZERO=1e-10;         // less precision than in Instance
    const double MY_ONE=1-MY_ZERO;
    int N=numNodes()-1;
    // depot has no demand
    probs.push_back(vector<double>());
    for (int i=1; i<=N; ++i) {
        vector<double> iprobs;
        double cum=0;
        int k=0;
        // adjustment for numerically unstable instances
        int adj=1;
        if (instance.id()==5003 || instance.id()==5004)
            adj=1e+1;
        while (cum<MY_ONE/adj) {
            double p=instance.prob(i, k);
            iprobs.push_back(p);
            cum+=p;
            k++;
        }
        probs.push_back(iprobs);
    }
    // TODO: make threshold work in a way that probs sum up to one
    // setup maxdemands and mindemands
    maxdemands.push_back(-1);
    mindemands.push_back(-1);
    for (int i=1; i<=N; ++i) {
        for (int k=0; k<probs[i].size(); ++k) {
            if (probs[i][k]>=Config::THRESH_PROB_ZERO) {
                mindemands.push_back(k);
                break;
            }
        }
        for (int k=probs[i].size()-1; k>=0; --k) {
            if (probs[i][k]>=Config::THRESH_PROB_ZERO) {
                maxdemands.push_back(k);
                if (k>=3*instance.vehicleQ()) {
                    cerr<<"failed assumption: max demand >= 3*Q"<<endl;
                    cerr<<"(correct gamma() optimizations)"<<endl;
                    exit(-1);
                }
                break;
            }
        }
    }
    // sanity check
    for (int i=1; i<=N; ++i) {
        cout<<"customer "<<i<<" demands in ["<<mindemands[i]<<","
                <<maxdemands[i]<<"]"<<endl;
        if (mindemands[i]>maxdemands[i]) {
            cerr<<"mindemands["<<i<<"]=="<<mindemands[i]
                    <<" and maxdemands["<<i<<"]=="<<maxdemands[i]<<endl;
            exit(-1);
        }
    }
}

/* Works for loads up to 5 */
void Data::genMinDists() {
    for (int i=0; i<numNodes(); ++i) {
        // min dist from node i (can be the depot)
        vector<double> dists;
        for (int l=0; l<=5*instance.vehicleQ(); ++l) {
            // ...to some node j such that expDemand(j)<=l
            double mdist=1024*1024;     // some large number
            for (int j=1; j<numNodes(); ++j)
                if (j!=i && expDemand(j)<=l && instance.dist(i, j)<mdist)
                    mdist=instance.dist(i, j);
            dists.push_back(mdist);
        }
        mindists.push_back(dists);
    }
}

void Data::genMinDistsPlus() {
    struct aux {
        int idx;
        double dist;
        aux(int i, double d) : idx{i}, dist{d} {}
    };
    for (int i=0; i<numNodes(); ++i) {
        vector<aux> vaux;
        for (int j=1; j<numNodes(); ++j) {
            if (j!=i) {
                vaux.push_back(aux(j, dist(i, j)));
            }
        }
        sort(vaux.begin(), vaux.end(),
                [](const aux& a1, const aux& a2) {return a1.dist<a2.dist;});
        vector<int> dists;
        for (auto& a : vaux)
            dists.push_back(a.idx);
        mindistsplus.push_back(dists);
    }
    cout<<"mindistsplus:"<<endl;
    for (int i=0; i<numNodes(); ++i) {
        cout<<i<<": ";
        for (auto j : mindistsplus.at(i))
            cout<<j<<" ";
        cout<<endl;
    }
}

/* Generate initial (maybe infeasible) routes so CG can start */
vector<Route> Data::initialRoutes() {
    vector<Route> routes;
    {
        Route giant {0};
        int load=0;
        for (int i=1; i<numNodes(); ++i) {
            giant.addNode(i);
            load+=expDemand(i);
        }
        giant.addNode(0);
        if (load<=maxLoad(true)) {
            cerr<<"giant tour is feasible! implement TSP initial route!"<<endl;
            exit(-1);
        }
        giant.setLoad(load);
        PolicySolverI sdp(instance, giant);
        double ecost=sdp.solve();
        cout<<"expected cost of giant tour: "<<ecost<<endl;
        routes.push_back(giant);
    }
    // and the depot-node-depot routes {0,i,0}
    for (int i=1; i<numNodes(); ++i) {
        Route dnd {0};
        dnd.addNode(i);
        dnd.addNode(0);
        if (expDemand(i)>maxLoad(true)) {
            cerr<<"expected demand of "<<i<<" is greater than maxLoad"<<endl;
            exit(-1);
        }
        dnd.setLoad(expDemand(i));
        PolicySolverI sdp(instance, dnd);
        double ecost=sdp.solve();
        dnd.setExpCost(ecost);
        routes.push_back(dnd);
        cout<<dnd<<endl;
    }
    return routes;
}

int Data::minVehicles(const vector<int>& nodes) const {
    int sumd=0;
    for (auto i : nodes)
        sumd+=expDemand(i);
    return std::ceil((sumd*1.0)/maxLoad(true));
}

vector<Route> Data::optCVRPRoutes() const {
    vector<Route> routes;
    vector<vector<int>> cvrpRoutes=instance.optCVRPRoutes();
    for (const auto& r : cvrpRoutes)
        routes.push_back(Route(r));
    return routes;
}

void Data::printConfigParams() const {
    cout<<"INITIAL_UB="<<Config::INITIAL_UB<<endl;
    cout<<"BRANCH_COEF_INT="<<Config::BRANCH_COEF_INT<<endl;
    cout<<"EPS_BINARY_SOL="<<Config::EPS_BINARY_SOL<<endl;
    cout<<"INFEASIBLE_ROUTE_PENALTY="<<Config::INFEASIBLE_ROUTE_PENALTY<<endl;
    cout<<"COLS_LIMIT_EXACT="<<Config::COLS_LIMIT_EXACT<<endl;
    cout<<"COLS_LIMIT_HEUR="<<Config::COLS_LIMIT_HEUR<<endl;
    cout<<"CG_ACCEL_RATIO_MAXLOAD="<<Config::CG_ACCEL_RATIO_MAXLOAD<<endl;
    cout<<"CG_ACCEL_ITERATIONS="<<Config::CG_ACCEL_ITERATIONS<<endl;
    cout<<"THREADS_LSOLVER="<<Config::THREADS_LSOLVER<<endl;
    cout<<"HEURISTIC="<<Config::HEURISTIC<<endl;
    cout<<"MAX_MEM_LABELS="<<Config::MAX_MEM_LABELS<<endl;
    cout<<"ENABLE_CUTS="<<Config::ENABLE_CUTS<<endl;
    cout<<"MAX_CUTS="<<Config::MAX_CUTS<<endl;
    cout<<"EPS_EDGE_NONZERO="<<Config::EPS_EDGE_NONZERO<<endl;
    cout<<"RCC_LEQ="<<Config::RCC_LEQ<<endl;
    cout<<"MIN_TO_PRINT="<<Config::MIN_TO_PRINT<<endl;
    cout<<"THRESH_PROB_ZERO="<<Config::THRESH_PROB_ZERO<<endl;
}

/* Profiling */
void Data::endProfilingAddNode() {
    chrono::steady_clock::time_point end=std::chrono::steady_clock::now();
    elapAddNode+=chrono::duration_cast<chrono::microseconds>
            (end-beginAddNode).count();
}

void Data::endProfilingAddNodeDirectCost() {
    chrono::steady_clock::time_point end=std::chrono::steady_clock::now();
    elapAddNodeDirectCost+=chrono::duration_cast<chrono::microseconds>
            (end-beginAddNodeDirectCost).count();
}

void Data::memoryProfilingLabels(long cap, long used) {
    double mem=sizeULabel*used+sizeRLabel*(cap-used);
    if (mem>maxMemoryLabels) {
        maxMemoryLabels=mem;
        cout<<"profiling: active labels: "<<used<<" ("<<(used*100.0)/cap
                <<"\%)\tmemory: "<<mem/(1024*1024*1024)<<"G"<<endl;
        if (mem/(1024*1024*1024)>Config::MAX_MEM_LABELS) {
            cout<<"swap imminent, terminating execution"<<endl;
            exit(-1);
        }
    }
}

void Data::printProfiling() {
    cout<<"TIME PROFILING"<<endl;
    cout<<"adding nodes: "<<elapAddNode/1e+6<<" secs"<<endl;
    cout<<"\t... dcost: "<<elapAddNodeDirectCost/1e+6<<" secs"<<endl;
    cout<<"MEMORY PROFILING"<<endl;
    cout<<"profiling: peak memory allocated to labels: "
            <<maxMemoryLabels/(1024*1024*1024)<<"G"<<endl;
}

void Data::startProfilingAddNode() {
    beginAddNode=chrono::steady_clock::now();
}

void Data::startProfilingAddNodeDirectCost() {
    beginAddNodeDirectCost=chrono::steady_clock::now();
}

