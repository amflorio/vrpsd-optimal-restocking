#include <vector>
#include "cvrpsep/basegrph.h"
#include "cvrpsep/cnstrmgr.h"
#include "cvrpsep/capsep.h"
#include "Config.h"
#include "SeparationRCC.h"

using namespace std;

shared_ptr<Data> SeparationRCC::data;

void SeparationRCC::checkCut(const RoundedCapacityCut& cut,
        const vector<int>& nodes) const {
    double sumaijs=0.0;
    for (const auto& a : cut.arcs)
        sumaijs+=a.val;
    // nodes is a binary vector, convert it to right format before call below
    vector<int> vnodes;
    for (int i=0; i<nodes.size(); ++i)
        if (nodes[i])
            vnodes.push_back(i);
    int minV=data->minVehicles(vnodes);
    if (Config::RCC_LEQ) {
        if (sumaijs<=vnodes.size()-minV) {
            cout<<"checkCut(): inconsistency: cut is not violated"<<endl;
            cout<<"checkCut(): |S|="<<vnodes.size()<<"\tsumaijs="<<sumaijs
                    <<"\tk(S)="<<minV<<endl;
            exit(-1);
        }
    } else {
        if (sumaijs>=minV) {
            cout<<"checkCut(): inconsistency: cut is not violated"<<endl;
            cout<<"checkCut(): |S|="<<vnodes.size()<<"\tsumaijs(delta(S))="
                    <<sumaijs<<" >= "<<minV<<"=k(S)"<<endl;
            exit(-1);
        }
    }
}

vector<RoundedCapacityCut> SeparationRCC::solve() const {
    int N=data->numNodes();     // # of customers: N-1
    // setting up demands
    vector<int> demands {0};
    for (int i=1; i<N; ++i)
        demands.push_back(data->expDemand(i));
    // setting up edge values (we need edges for the cvrpsep)
    vector<vector<double>> xijs(N, vector<double>(N, 0.0));
    vector<vector<double>> aijs(N, vector<double>(N, 0.0));
    for (const auto& rs : sol.routesSolution()) {
        vector<int> nodes=rs.route.nodes();
        for (int k=1; k<nodes.size(); ++k) {
            // i<=j
            aijs[nodes[k-1]][nodes[k]]+=rs.coef;
            int i=nodes[k-1]<nodes[k]?nodes[k-1]:nodes[k];
            int j=nodes[k]>nodes[k-1]?nodes[k]:nodes[k-1];
            xijs[i][j]+=rs.coef;
        }
    }
    // setting up the necessary edge data structures for cvrpsep call
    vector<int> edgex {0};
    vector<int> edgey {0};          // edge numbers are 1-based (see manual)
    vector<double> edgeval {0.0};
    for (int i=0; i<N; ++i)
        for (int j=i+1; j<N; ++j)
            if (xijs[i][j]>Config::EPS_EDGE_NONZERO) {
                edgex.push_back(i==0?N:i);      // see cvrpsep manual
                edgey.push_back(j);
                edgeval.push_back(xijs[i][j]);
            }
    CnstrMgrPointer cutsCMP=nullptr, oldCutsCMP=nullptr;
    CMGR_CreateCMgr(&cutsCMP, 100);
    CMGR_CreateCMgr(&oldCutsCMP, 100);
    char intAndFeasible;
    double maxViolation;
    CAPSEP_SeparateCapCuts(N-1, demands.data(), data->maxLoad(true),
            edgex.size()-1, edgex.data(), edgey.data(), edgeval.data(),
            oldCutsCMP, Config::MAX_CUTS, 1e-4, &intAndFeasible, &maxViolation,
            cutsCMP);
    //cout<<cutsCMP->Size<<" cut(s) found!"<<endl;
    if (cutsCMP->Size>0)
        cout<<"max violation: "<<maxViolation<<endl;
    /*
    if (intAndFeasible)
        cout<<"integer and feasible"<<endl;
    */
    // retrieving cuts
    vector<RoundedCapacityCut> cuts;
    for (int i=0; i<cutsCMP->Size; ++i) {
        vector<int> cutnodes(data->numNodes(), 0);
        for (int j=1; j<=cutsCMP->CPL[i]->IntListSize; ++j) {
            int n=cutsCMP->CPL[i]->IntList[j];
            cutnodes.at(n)=1;
            // assumption: depot is never included in any cut
            if (n==0 || n==N) {     // manual not clear
                cout<<"assertion failed: depot is included in cut"<<endl;
                /*
                cout<<"if some arc to/from the depot is included in a cut, "<<
                        "pricing must be changed accordingly (addRootLabels, "<<
                        "possibly others), and checking for a route with "<<
                        "negative red. cost (currently in Label) must take "<<
                        "into account eventual arcs to/from the depot"<<endl;
                */
                exit(-1);
            }
        }
        assert(cutnodes[0]==0);
        vector<Arc> arcs;
        if (Config::RCC_LEQ) {
            for (int ai=0; ai<cutnodes.size(); ++ai) {
                assert(aijs[ai][ai]==0);
                if (cutnodes[ai])
                    for (int aj=0; aj<cutnodes.size(); ++aj)
                        if (cutnodes[aj] && aj!=ai)
                            arcs.push_back(Arc(ai, aj, aijs[ai][aj]));
            }
            // CPL[i]->RHS gives |S|-k(S)
            int rhs=cutsCMP->CPL[i]->RHS;
            cuts.push_back(RoundedCapacityCut(arcs, rhs));
        } else {
            for (int ai=0; ai<cutnodes.size(); ++ai)
                if (cutnodes[ai])
                    for (int aj=0; aj<cutnodes.size(); ++aj)
                        if (!cutnodes[aj])
                            arcs.push_back(Arc(ai, aj, aijs[ai][aj]));
            // CPL[i]->RHS gives |S|-k(S), but we need only k(S)
            int rhs=cutsCMP->CPL[i]->IntListSize-cutsCMP->CPL[i]->RHS;
            cuts.push_back(RoundedCapacityCut(arcs, rhs));
        }
        checkCut(cuts.back(), cutnodes);
    }
    return cuts;
}

