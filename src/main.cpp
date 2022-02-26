#include <chrono>
#include <iomanip>
#include <iostream>
#include "Config.h"
#include "Instance.h"
#include "IntegerSolver.h"
#include "Label.h"
#include "LinearSolver.h"
#include "PolicySolverI.h"
#include "Pricing.h"
#include "Route.h"
#include "Separation3SRC.h"
#include "SeparationRCC.h"
#include "Solution.h"
#include "Solver.h"
#include "SwitchPolicy.h"

using namespace std;

void mycontinue() {
    cout<<"press <Enter> to continue"<<endl;
    cin.get();
}

int main(int argc, char* argv[]) {
    chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    if (argc!=3) {
        cerr<<"usage: "<<argv[0]<<" <mode> <instance ID>"<<endl;
        cerr<<"where <mode> ="<<endl;
        cerr<<"\t0   Branch and Price"<<endl;
        cerr<<"\t1   Comparison Opt. Restocking and Switch Policies"<<endl;
        return 0;
    }
    int mode=std::stoi(argv[1]);
    int id=std::stoi(argv[2]);
    cout<<fixed;
    cout<<setprecision(3);
    Instance inst(id);
    inst.printInfo();
    if (mode==0) {
        cout<<"*** BRANCH AND PRICE ***"<<endl;
        shared_ptr<Data> data {new Data(inst)};
        // share the common problem data with all the BnP modules that need it
        Label::setData(data);
        Pricing::setData(data);
        Route::setData(data);
        SeparationRCC::setData(data);
        Separation3SRC::setData(data);
        SubsetRowCut::setData(data);
        IntegerSolver::setData(data);
        LinearSolver::setData(data);
        Solution::setData(data);
        Solver::setData(data);
        // deterministic equivalent
        cout<<"Deterministic equivalent solution:"<<endl;
        double detEq=0.0;
        double apriori=0.0;
        vector<Route> cvrpOpt=data->optCVRPRoutes();
        for (Route& one : cvrpOpt) {
            apriori+=one.aPriori();
            // we analyze both orientations of each optimal CVRP route
            double onecost;
            {
                PolicySolverI sdp(inst, one);
                onecost=sdp.solve();
                one.setExpCost(onecost);
            }
            Route rev=one.reverse();
            double revcost;
            {
                PolicySolverI sdp(inst, rev);
                revcost=sdp.solve();
                rev.setExpCost(revcost);
            }
            cout<<(onecost<revcost?one:rev)<<endl;
            detEq+=(onecost<revcost?onecost:revcost);
        }
        cout<<"a-priori length of the DE solution: "<<apriori<<endl;
        cout<<"deterministic equivalent solution value: "<<detEq<<endl;
        // initial upper-bound for the branch-and-bound
        Solution ub;
        ub.setValue(Config::INITIAL_UB);
        Solver solver(vector<Column>(), ub);
        Solution sol=solver.solve();
        data->printConfigParams();
        sol.print();
        data->printStats();
        cout<<"a-priori length of the DE solution: "<<apriori<<endl;
        cout<<"Value of the Stochastic Solution (VSS) = "<<(detEq-sol.value())
                <<"\t("<<100*((detEq-sol.value())/detEq)<<"\%)"<<endl;
        cout<<"BnP solution (value): "<<sol.value()<<endl;
        cout<<"BnP solution (routes):"<<endl;
        for (auto& r : sol.routes())
            cout<<r<<endl;
        cout<<"verifying route expected costs with the SDP algorithm"<<endl;
        SwitchPolicy::setData(data);
        double totOR=0.0;
        double totSP=0.0;
        for (auto& r : sol.routes()) {
            PolicySolverI sdp(inst, r);
            double sdpcost=sdp.solve();
            totOR+=sdpcost;
            cout<<"Route: "<<r.expCost()<<"\tSDP: "<<sdpcost<<"\toff by: "
                    <<((r.expCost()-sdpcost)/sdpcost)*100<<"\%"<<endl;
            SwitchPolicy sp(r);
            double sc=sp.solve();
            totSP+=sc;
            double diff=((sc-sdpcost)/sdpcost)*100;
            cout<<"(under Switch policy: "<<sc<<"\t"<<diff<<"\%)"<<endl;
        }
        double diff=((totSP-totOR)/totOR)*100;
        cout<<"Switch policy (total): "<<diff<<"\%"<<endl;
        data->printProfiling();
    } else if (mode==1) {
        // redirect "normal" output so the export output comes clean
        streambuf *old=cout.rdbuf();
        stringstream ss;
        cout.rdbuf(ss.rdbuf());
        shared_ptr<Data> data {new Data(inst)};
        SwitchPolicy::setData(data);
        Route::setData(data);
        vector<Route> cvrpRoutes=data->optCVRPRoutes();
        cout.rdbuf(old);
        double totOR=0.0;
        double totSP=0.0;
        for (const auto& r : cvrpRoutes) {
            cout<<r<<endl;
            PolicySolverI orp(inst, r);
            double orc=orp.solve();
            Route rev=r.reverse();
            PolicySolverI orprev(inst, rev);
            double orcrev=orprev.solve();
            if (orcrev<orc)
                orc=orcrev;
            totOR+=orc;
            cout<<"Opt restocking: "<<orc<<endl;
            SwitchPolicy sp(r);
            double sc=sp.solve();
            SwitchPolicy sprev(rev);
            double screv=sprev.solve();
            if (screv<sc)
                sc=screv;
            totSP+=sc;
            double diff=((sc-orc)/orc)*100;
            cout<<"Switch policy: "<<sc<<"\t\t("<<diff<<"\%)"<<endl;
        }
        double diff=((totSP-totOR)/totOR)*100;
        cout<<"Switch policy (total): "<<diff<<"\%"<<endl;
    } else {
        cerr<<"invalid mode"<<endl;
        return 1;
    }
    chrono::steady_clock::time_point end=std::chrono::steady_clock::now();
    int elap=chrono::duration_cast<chrono::seconds>(end-begin).count();
    cout<<"time elapsed: "<<elap<<" secs"<<endl;
    return 0;
}

