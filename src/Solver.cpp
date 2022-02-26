#include <algorithm>
#include <cmath>
#include <memory>
#include "Config.h"
#include "IntegerSolver.h"
#include "LinearSolver.h"
#include "Route.h"
#include "Solution.h"
#include "Solver.h"

using namespace std;

shared_ptr<Data> Solver::data;

void Solver::addNewColumnsCuts(const vector<Column>& newcols) {
    for (const auto& col : newcols) {
        assert(col.coefsRCC.size()==RCCs.size());
        assert(col.coefsSRC.size()==SRCs.size());
        assert(col.compliesWithBranching);
        columns.push_back(col);
    }
    // update cut coefficients of all columns that need updating
    for (auto& col : columns) {
        // RCCs
        assert(col.coefsRCC.size()<=RCCs.size());
        for (int i=col.coefsRCC.size(); i<RCCs.size(); ++i)
            col.addRCCCoef(RCCs[i]);
        // SRCs
        assert(col.coefsSRC.size()<=SRCs.size());
        for (int i=col.coefsSRC.size(); i<SRCs.size(); ++i)
            col.addSRCCoef(SRCs[i]);
    }
}

bool Solver::assertCandidatesNotBranched(const vector<Arc>& candidates,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) const {
    //cout<<"assertCandidatesNotBranched..."<<endl;
    for (auto& c : candidates) {
        for (auto& a : arcsOn) {
            if (a.i==c.i && a.j==c.j) {
                cerr<<"assertion failed: candidate ["<<c.i<<","<<c.j
                        <<"] already branched on"<<endl;
                return false;
            }
        }
        for (auto& a : arcsOff) {
            if (a.i==c.i && a.j==c.j) {
                cerr<<"assertion failed: candidate ["<<c.i<<","<<c.j
                        <<"] already branched off"<<endl;
                return false;
            }
        }
    }
    return true;
}

void Solver::assertSolutionValid(const Solution& sol, const vector<Arc>& arcsOn,
        const vector<Arc>& arcsOff) const {
    cout<<"asserting solution is valid... "<<arcsOn.size()<<" arcs on , "
            <<arcsOff.size()<<" arcs off"<<endl;
    for (auto& rs : sol.routesSolution()) {
        if (rs.coef<=1e-5)
            continue;
        for (auto& arc : arcsOn) {
            if (!rs.route.complyWithArcOn(arc)) {
                cerr<<"assertion failed: arc "<<arc<<" is on, but the "<<
                        "following route is not compliant:"<<endl;
                cerr<<rs.coef<<" "<<rs.route<<endl;
                cerr<<"solution, arcs on and arcs off:"<<endl;
                sol.print();
                printArcsOn(arcsOn);
                printArcsOff(arcsOff);
                exit(-1);
            }
            assert(rs.route.complyWithArcOn(arc));
        }
        for (auto& arc : arcsOff) {
            assert(rs.route.complyWithArcOff(arc));
            /*
            if (!rs.route.complyWithArcOff(arc)) {
                cerr<<"assertion failed: arc "<<arc<<" is off, but the "<<
                        "following route is not compliant:"<<endl;
                cerr<<rs.coef<<" "<<rs.route<<endl;
                cerr<<"solution, arcs on and arcs off:"<<endl;
                sol.print();
                printArcsOn(arcsOn);
                printArcsOff(arcsOff);
                exit(-1);
            }
            */
        }
    }
}

/* Asserts the arc being branched on is valid. A valid arc has not already been
 * branched on. */
void Solver::assertValidArcToBranch(const Arc& arc, const vector<Arc>& arcsOn,
        const vector<Arc>& arcsOff) const {
    for (auto& a : arcsOn)
        assert(a.i!=arc.i || a.j!=arc.j);
    for (auto& a : arcsOff)
        assert(a.i!=arc.i || a.j!=arc.j);
}

void Solver::checkUpperBoundImproved(const LinearSolver& lsolver,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) {
    if (lsolver.hasIntermediateBinary()) {
        Solution interm=lsolver.getIntermediateBinary();
        if (interm.value()<upperBound.value()) {
            //assertSolutionValid(interm, arcsOn, arcsOff);
            updateUpperBound(interm);
        }
    }
}

/* Returns a vector of Arc where branching can be performed. The vector is
 * guaranteed to contain at least the Arc with fractional value closest to
 * 0.5. */
vector<Arc> Solver::arcToBranchCandidates(const Solution& sol) {
    vector<vector<double>> sumcoeffs(data->numNodes(),
            vector<double>(data->numNodes(), 0.0));
    // TODO: replace below by the more efficient filling loop in Separation
    for (auto& rs : sol.routesSolution())
        for (int i=1; i<data->numNodes(); ++i)
            for (int j=1; j<data->numNodes(); ++j)
                if (i!=j && rs.route.timesArc(i, j))
                    sumcoeffs[i][j]+=rs.coef*rs.route.timesArc(i, j);
    vector<Arc> candidates;
    double fracval=-1;
    int closest_i=-1;
    int closest_j=-1;
    for (int i=1; i<data->numNodes(); ++i) {
        for (int j=1; j<data->numNodes(); ++j) {
            if (abs(sumcoeffs[i][j]-0.5)<abs(fracval-0.5)) {
                fracval=sumcoeffs[i][j];
                closest_i=i;
                closest_j=j;
            }
            if (sumcoeffs[i][j]>=Config::BRANCH_COEF_INT
                    && sumcoeffs[i][j]<=1-Config::BRANCH_COEF_INT)
                candidates.push_back(Arc(i, j, sumcoeffs[i][j]));
        }
    }
    // arc with sumcoeffs closest to 0.5 is always a candidate
    if (candidates.size()==0)
        candidates.push_back(Arc(closest_i, closest_j, fracval));
    return candidates;
}

vector<Column> Solver::filterColsBranching(const vector<Column>& cols,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) const {
    cout<<"filtering columns for branching: cols.size()="<<cols.size()<<endl;
    vector<Column> filtered=cols;
    for (auto& arc : arcsOn)
        filtered=filterColsArcOn(filtered, arc);
    for (auto& arc : arcsOff)
        filtered=filterColsArcOff(filtered, arc);
    cout<<"filtering columns: filtered.size()="<<filtered.size()<<endl;
    return filtered;
}

vector<Column> Solver::filterCols(vector<Column> cols,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) const {
    for (auto& c : cols) {
        assert(c.compliesWithBranching);
        for (auto& arc : arcsOn)
            if (!c.route.complyWithArcOn(arc))
                c.compliesWithBranching=false;
        for (auto& arc : arcsOff)
            if (!c.route.complyWithArcOff(arc))
                c.compliesWithBranching=false;
    }
    // I now this can be more efficient but this is not relevant
    return cols;
}

vector<int> Solver::filterColsIndex(const vector<Column>& cols,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) const {
    vector<int> filtered(cols.size(), 1);
    for (int i=0; i<cols.size(); ++i) {
        if (cols[i].expCost==Config::INFEASIBLE_ROUTE_PENALTY)
            continue;
        for (auto& a : arcsOn) {
            if (!cols[i].route.complyWithArcOn(a)) {
                filtered[i]=0;
                break;
            }
        }
        if (filtered[i]==1) {
            for (auto& a : arcsOff) {
                if (!cols[i].route.complyWithArcOff(a)) {
                    filtered[i]=0;
                    break;
                }
            }
        }
    }
    return filtered;
}

vector<Column> Solver::filterColsArcOn(const vector<Column>& cols,
        const Arc& arc) const {
    vector<Column> filtered;
    for (auto& c : cols) {
        // always keep the initial (infeasible) route(s)
        if (c.expCost==Config::INFEASIBLE_ROUTE_PENALTY)
            filtered.push_back(c);
        else if (c.route.complyWithArcOn(arc))
            filtered.push_back(c);
    }
    return filtered;
}

vector<Column> Solver::filterColsArcOff(const vector<Column>& cols,
        const Arc& arc) const {
    vector<Column> filtered;
    for (auto& c : cols) {
        // always keep the initial (infeasible) route(s)
        if (c.expCost==Config::INFEASIBLE_ROUTE_PENALTY)
            filtered.push_back(c);
        else if (c.route.complyWithArcOff(arc))
            filtered.push_back(c);
    }
    return filtered;
}

BnBNode Solver::nextNode() {
    double minval=1024.0*1024.0;
    int minidx=-1;
    for (int i=0; i<pending.size(); ++i) {
        if (pending[i].sol.value()<minval) {
            minval=pending[i].sol.value();
            minidx=i;
        }
    }
    BnBNode next=pending[minidx];
    pending.erase(pending.begin()+minidx);
    if (!Config::HEURISTIC && minval+1e-6<lowerBound) {
        cerr<<setprecision(6);
        cout<<setprecision(6);
        cerr<<"inconsistency: LB dec from "<<lowerBound<<" to "<<minval<<endl;
        cerr<<"diff*1e+6="<<(lowerBound-minval)*1e+6<<endl;
        exit(-1);
    }
    lowerBound=minval;
    return next;
}

void Solver::printArcsOn(const vector<Arc>& arcsOn) const {
    cout<<"arcs on:";
    for (auto& arc : arcsOn)
        cout<<" "<<arc;
    cout<<endl;
}

void Solver::printArcsOff(const vector<Arc>& arcsOff) const {
    cout<<"arcs off:";
    for (auto& arc : arcsOff)
        cout<<" "<<arc;
    cout<<endl;
}

/* Returns the Arc which branches to the highest (heuristic) lower-bound. */
Arc Solver::selectBestCandidate(const vector<Arc>& candidates,
        const vector<Arc>& arcsOn, const vector<Arc>& arcsOff) {
    double highest=0;
    int high_i=-1;
    int high_j=-1;
    double high_val=-1;
    vector<Column> cols=filterColsBranching(columns, arcsOn, arcsOff);
    // solveBranching does not invoke Pricing, so no need for arcsOn or arcsOff
    LinearSolver lsolver(cols, RCCs, SRCs, vector<Arc>(), vector<Arc>());
    for (auto& arc : candidates) {
        // with soft-branching, an arc could be chosen more than once
        bool skip=false;
        for (auto& a : arcsOn)
            if (a.i==arc.i && a.j==arc.j)
                skip=true;
        for (auto& a : arcsOff)
            if (a.i==arc.i && a.j==arc.j)
                skip=true;
        if (skip)
            continue;       // make sure an arc is not chosen more than once
        double lbon;
        double lboff;
        {
            vector<Arc> newArcsOn;
            newArcsOn.push_back(arc);
            Solution sol=lsolver.solveBranching(filterColsIndex(cols, newArcsOn,
                    vector<Arc>()));
            lbon=sol.value();
        }
        {
            vector<Arc> newArcsOff;
            newArcsOff.push_back(arc);
            Solution sol=lsolver.solveBranching(filterColsIndex(cols,
                    vector<Arc>(), newArcsOff));
            lboff=sol.value();
        }
        double lbmin=min(lbon, lboff);
        cout<<"lb (heur) inc to "<<lbmin<<" when branching on "<<arc
                <<" with a frac val "<<arc.val<<endl;
        if (lbmin>highest) {
            high_i=arc.i;
            high_j=arc.j;
            high_val=arc.val;
            highest=lbmin;
        }
    }
    assert(high_i!=-1);
    return Arc(high_i, high_j, high_val);
}

Solution Solver::solve() {
    auto isolver=make_shared<IntegerSolver>();
    {
        // root node
        chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
        LinearSolver lsolver(columns, RCCs, SRCs, vector<Arc>(), vector<Arc>());
        Solution sol=lsolver.solveHeuristic();
        sol=lsolver.solveExact(true, upperBound.value(), isolver);
        if (!sol.feasible()) {
            cout<<"bnb: instance "<<data->getInstanceID()<<endl;
            cout<<"bnb: root node infeasible"<<endl;
            cout<<"bnb: maxLoad = "<<data->maxLoad(true)<<endl;
            exit(-1);
        }
        lowerBound=sol.value();
        cout<<"bnb: root node: lower bound = "<<lowerBound<<endl;
        chrono::steady_clock::time_point end=std::chrono::steady_clock::now();
        int elap=chrono::duration_cast<chrono::seconds>(end-begin).count();
        cout<<"bnb: root node: elapsed time: "<<elap<<" secs"<<endl;
        data->printProfiling();
        //exit(-1);
        checkUpperBoundImproved(lsolver, vector<Arc>(), vector<Arc>());
        vector<Column> cols=lsolver.getNewColumns();
        /*
        if (Config::FIX_VEHICLES) {
            cols.erase(std::remove_if(cols.begin(), cols.end(),
                    [] (const Column& c) {return !c.onceUsed;}), cols.end());
        }
        */
        addNewColumnsCuts(cols);
        //addNewColumnsCuts(lsolver.getNewColumns());
        pending.push_back(BnBNode(sol, vector<Arc>(), vector<Arc>()));
    }
    while (pending.size()>0) {
        Solution isol=isolver->solve(columns);
        if (isol.feasible() && isol.value()<upperBound.value()) {
            cout<<"bnb: new upper-bound obtained by the int solver"<<endl;
            updateUpperBound(isol);
        }
        cout<<"bnb: "<<pending.size()<<" pending nodes"<<endl;
        BnBNode n=nextNode();
        cout<<"bnb: selected node with value = "<<n.sol.value()
                <<"   upperBound.value() = "<<upperBound.value()<<endl;
        cout<<"bnb: optimality gap: "
                <<100*((upperBound.value()-lowerBound)/lowerBound)<<"\%"<<endl;
        if (upperBound.value()-lowerBound<1e-3)
            return upperBound;
        if (!n.sol.feasible()) {
            cerr<<"bnb: inconsistency: infeasible node"<<endl;
            exit(-1);
        }
        if (n.sol.value()>upperBound.value()) {
            cout<<"bnb: pruning node: upper bound unachievable"<<endl;
            continue;
        }
        if (n.sol.binary()) {
            if (n.sol.value()+1e-5<upperBound.value()) {
                cerr<<"bnb: inconsistency: binary sol. with better value"<<endl;
                exit(-1);
            }
            cout<<"bnb: pruning node: upper bound unachievable"<<endl;
            continue;
        }
        // branch
        vector<Arc> candidates=arcToBranchCandidates(n.sol);
        /*
        if (!assertCandidatesNotBranched(candidates, n.arcsOn, n.arcsOff)) {
            cerr<<"solution:"<<endl;
            n.sol.print();
            printArcsOn(n.arcsOn);
            printArcsOff(n.arcsOff);
            exit(-1);
        }
        */
        Arc arc=selectBestCandidate(candidates, n.arcsOn, n.arcsOff);
        assertValidArcToBranch(arc, n.arcsOn, n.arcsOff);
        cout<<"branching on arc "<<arc<<", frac val of "<<arc.val<<endl;
        {
            vector<Arc> newArcsOn {n.arcsOn};
            newArcsOn.push_back(arc);
            LinearSolver lsolver(filterCols(columns, newArcsOn, n.arcsOff),
                    RCCs, SRCs, newArcsOn, n.arcsOff);
            lsolver.solvePRA();
            Solution sol=lsolver.solveHeuristic();
            sol=lsolver.solveExact(true, upperBound.value(), nullptr);
            //assertSolutionValid(sol, newArcsOn, n.arcsOff);
            addNewColumnsCuts(lsolver.getNewColumns());
            /*
            if (sol.feasible())
                assertSolutionValid(sol, newArcsOn, n.arcsOff);
            */
            if (sol.value()>upperBound.value())
                cout<<"bnb: pruning node when branching"<<endl;
            else {
                checkUpperBoundImproved(lsolver, newArcsOn, n.arcsOff);
                if (!sol.binary())
                    pending.push_back(BnBNode(sol, newArcsOn, n.arcsOff));
                else if (sol.value()+1e-5<upperBound.value()) {
                    cerr<<"bnb: inconsistency: binary sol. not updated"<<endl;
                    exit(-1);
                }
            }
        }
        {
            vector<Arc> newArcsOff {n.arcsOff};
            newArcsOff.push_back(arc);
            LinearSolver lsolver(filterCols(columns, n.arcsOn, newArcsOff),
                    RCCs, SRCs, n.arcsOn, newArcsOff);
            lsolver.solvePRA();
            Solution sol=lsolver.solveHeuristic();
            sol=lsolver.solveExact(true, upperBound.value(), nullptr);
            //assertSolutionValid(sol, n.arcsOn, newArcsOff);
            addNewColumnsCuts(lsolver.getNewColumns());
            /*
            if (sol.feasible())
                assertSolutionValid(sol, n.arcsOn, newArcsOff);
            */
            if (sol.value()>upperBound.value())
                cout<<"bnb: pruning node when branching"<<endl;
            else {
                checkUpperBoundImproved(lsolver, n.arcsOn, newArcsOff);
                if (!sol.binary())
                    pending.push_back(BnBNode(sol, n.arcsOn, newArcsOff));
                else if (sol.value()+1e-5<upperBound.value()) {
                    cerr<<"bnb: inconsistency: binary sol. not updated"<<endl;
                    exit(-1);
                }
            }
        }
    }
    cout<<"bnb: optimality gap: "
            <<100*((upperBound.value()-lowerBound)/lowerBound)<<"\%"<<endl;
    return upperBound;
}

void Solver::updateUpperBound(const Solution& newUB) {
    upperBound=newUB;
    cout<<"upper-bound improved: "<<upperBound.value()<<" (paths: "
            <<upperBound.routes().size()<<")"<<endl;
    upperBound.print();
}

