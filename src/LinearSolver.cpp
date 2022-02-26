#include <algorithm>
#include <exception>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include "Config.h"
#include "LinearSolver.h"
#include "Pricing.h"
#include "Separation3SRC.h"
#include "SeparationRCC.h"
#include "Solution.h"

using namespace std;

shared_ptr<Data> LinearSolver::data;

ILOSTLBEGIN

LinearSolver::LinearSolver(vector<Column> cols,
        vector<RoundedCapacityCut>& vRCCs, vector<SubsetRowCut>& vSRCs,
        vector<Arc> aOn, vector<Arc> aOff) : model{env}, cplex{env}, var{env},
        con{env}, obj{IloMinimize(env)}, columns{move(cols)}, RCCs{vRCCs},
        SRCs{vSRCs}, arcsOn{move(aOn)}, arcsOff{move(aOff)} {
    // adding RHS of max vehicles constraint
    if (Config::FIX_VEHICLES) {
        con.add(IloRange(env, -IloInfinity, data->maxVehicles()));
        typeConstraint.push_back(CTR_MAXV);
        // add surplus var to penalize violations of the max vehicles ctr
        IloNumVar surplus(env, 0, IloInfinity);
        con[0].setLinearCoef(surplus, -1);
        obj.setLinearCoef(surplus, 1024*1024);
    }
    // adding RHS of partitioning constraints
    for (int i=0; i<data->numNodes()-1; ++i) {
        con.add(IloRange(env, 1, 1));
        typeConstraint.push_back(CTR_PTT);
    }
    // adding RHS of RCCs
    for (const auto& c : RCCs) {
        con.add(IloRange(env, Config::RCC_LEQ?-IloInfinity:c.rhs,
                Config::RCC_LEQ?c.rhs:IloInfinity));
        //con.add(IloRange(env, c.rhs, IloInfinity));
        typeConstraint.push_back(CUT_RCC);
    }
    // adding RHS of SRCs
    for (const auto& c : SRCs) {
        con.add(IloRange(env, -IloInfinity, c.rhs()));
        typeConstraint.push_back(CUT_SRC);
    }
    model.add(obj);
    model.add(con);
    cplex.setParam(IloCplex::Threads, Config::THREADS_LSOLVER);
    cplex.setOut(env.getNullStream());
    cplex.extract(model);
    // the loaded columns will have to be included in the RMP
    for (auto& c : columns)
        c.included=false;
}

void LinearSolver::addRCCut(const RoundedCapacityCut& cut) {
    IloNumVarArray vars(env);
    IloNumArray vals(env);
    for (int i=0; i<columns.size(); ++i) {
        assert(columns[i].coefsRCC.size()==RCCs.size());
        vars.add(var[i]);
        vals.add(columns[i].addRCCCoef(cut));
    }
    IloRange ilocut(env, Config::RCC_LEQ?-IloInfinity:cut.rhs,
            Config::RCC_LEQ?cut.rhs:IloInfinity);
    //IloRange ilocut(env, cut.rhs, IloInfinity);
    ilocut.setLinearCoefs(vars, vals);
    con.add(ilocut);
    typeConstraint.push_back(CUT_RCC);
    model.add(ilocut);
    vars.end();
    vals.end();
    RCCs.push_back(cut);
}

void LinearSolver::addSRCut(const SubsetRowCut& cut) {
    IloNumVarArray vars(env);
    IloNumArray vals(env);
    for (int i=0; i<columns.size(); ++i) {
        assert(columns[i].coefsSRC.size()==SRCs.size());
        vars.add(var[i]);
        vals.add(columns[i].addSRCCoef(cut));
    }
    IloRange ilocut(env, -IloInfinity, cut.rhs());
    ilocut.setLinearCoefs(vars, vals);
    con.add(ilocut);
    typeConstraint.push_back(CUT_SRC);
    model.add(ilocut);
    vars.end();
    vals.end();
    SRCs.push_back(cut);
}

void LinearSolver::addNewColumn(Route route) {
    Column col(route, route.expCost(), false, true);
    // set onceUsed flag to the 0-i-0 routes (so they are always kept)
    if (route.numVisits()==1)
        col.onceUsed=true;
    // fill the coefficients of the RCCs
    for (const auto& c : RCCs)
        col.addRCCCoef(c);
    // and of the SRCs
    for (const auto& c : SRCs)
        col.addSRCCoef(c);
    columns.push_back(col);
}

vector<pair<SubsetRowCut, double>> LinearSolver::findNBestCuts(
        const vector<SubsetRowCut>& cuts, int N) const {
    vector<SubsetRowCut> newSRCs=SRCs;      // copy pre-existent SRCs
    newSRCs.insert(newSRCs.end(), cuts.begin(), cuts.end());
    LinearSolver ls(columns, RCCs, newSRCs, arcsOn, arcsOff);
    for (auto& col : ls.columns)
        for (const auto& c : cuts)
            col.addSRCCoef(c);
    const int os=Config::FIX_VEHICLES*1+data->numNodes()-1+RCCs.size()
            +SRCs.size();
    ls.solveRMP(false);     // so that the bounds set below 'stick'
    for (int k=0; k<cuts.size(); ++k) {
        assert(ls.typeConstraint[k+os]==CUT_SRC);
        ls.con[k+os].setBounds(-IloInfinity, IloInfinity);
    }
    vector<pair<SubsetRowCut, double>> srcs;
    for (int k=0; k<cuts.size(); ++k) {
        ls.con[k+os].setBounds(-IloInfinity, 1);    // TODO: only for 3-SRCs
        double val=ls.solveRMP(false).value();
        ls.con[k+os].setBounds(-IloInfinity, IloInfinity);
        cout<<"SRC: "<<cuts[k]<<" -> "<<val<<endl;
        srcs.emplace_back(cuts[k], val);
    }
    sort(srcs.begin(), srcs.end(), [](const pair<SubsetRowCut, double>& c1,
            const pair<SubsetRowCut, double>& c2){return c1.second>c2.second;});
    if (srcs.size()>N)
        srcs.erase(srcs.begin()+N, srcs.end());
    return srcs;
}

void LinearSolver::genInitialColumns() {
    for (auto& r : data->initialRoutes())
        addNewColumn(r);
    // generate a good set of columns with the acceleration algorithm
    int maxLoad=data->maxLoad(true);
    int l=maxLoad*Config::CG_ACCEL_RATIO_MAXLOAD+0.5;
    int inc=max(1, (maxLoad-l)/Config::CG_ACCEL_ITERATIONS);
    while (l<maxLoad) {
        data->setACGMaxLoad(l);
        Solution sol=solveAccelerated();
        cout<<"accel. CG: done for l="<<l<<"\tsol. value="<<sol.value()<<endl;
        l+=inc;
    }
    data->setACGMaxLoad(data->maxLoad(true));
    cout<<"ACG finished"<<endl;
    data->statSetACGFinished();
}

vector<Column> LinearSolver::getNewColumns() {
    vector<Column> newcols;
    for (auto& c : columns)
        if (c.isNew) {
            c.isNew=false;
            newcols.push_back(c);
        }
    return newcols;
}

/* Solves the linear relaxation of the VRPSD by column generation. This is a
 * simplified version of solveHeuristic, intended to be used only by the
 * accelerated column generation algorithm. */
Solution LinearSolver::solveAccelerated() {
    while (true) {
        Solution sol=solveRMP();
        Pricing pricing(sol, arcsOn, arcsOff, RCCs, SRCs);
        vector<Route> routes=pricing.solveHeuristic();
        for (auto& r : routes) {
            cout<<r<<endl;
            addNewColumn(r);
        }
        if (routes.size()<Config::COLS_LIMIT_HEUR)
            return solveRMP();
    }
}

/* Solves the linear relaxation of the VRPSD by column generation. Improves the
 * current solution until it becomes fractional, or until no further
 * improvement is possible. In this case an optimal binary solution is
 * returned. */
Solution LinearSolver::solveBinary() {
    data->statIncPricingSolveBinary();
    while (true) {
        currentSol=solveRMP();
        if (!currentSol.binary()) {
            cout<<"b: current solution is not binary"<<endl;
            return currentSol;
        }
        Pricing pricing(currentSol, arcsOn, arcsOff, RCCs, SRCs);
        vector<Route> routes=pricing.solveExact();
        if (routes.size()==0) {
            cout<<"b: no cols with neg rc found"<<endl;
            return currentSol;
        }
        cout<<"b: found "<<routes.size()<<" cols with neg rc:"<<endl;
        for (auto& r : routes) {
            cout<<r<<endl;
            addNewColumn(r);
        }
    }
}

Solution LinearSolver::solveBranching(const vector<int>& colidxs) {
    try {
        if (solveBranchingFirst) {
            for (auto& col : columns) {
                assert(col.compliesWithBranching);
                assert(!col.included);
                IloNumVar v(env, 0, IloInfinity);
                var.add(v);
                obj.setLinearCoef(v, col.expCost);
                if (Config::FIX_VEHICLES)
                    con[0].setLinearCoef(v, 1);
                for (int i=1; i<data->numNodes(); ++i)
                    con[Config::FIX_VEHICLES*1+i-1].setLinearCoef(v,
                            col.route.visitsNode(i));
                assert(col.coefsRCC.size()==RCCs.size());
                assert(col.coefsSRC.size()==SRCs.size());
                int offset=Config::FIX_VEHICLES*1+data->numNodes()-1;
                int next_rcc=0;
                int next_src=0;
                for (int i=0; i<RCCs.size()+SRCs.size(); ++i) {
                    assert(typeConstraint.at(i+offset)==CUT_RCC ||
                            typeConstraint[i+offset]==CUT_SRC);
                    if (typeConstraint[i+offset]==CUT_RCC) {
                        con[i+offset].setLinearCoef(v, col.coefsRCC[next_rcc]);
                        next_rcc++;
                    } else {
                        con[i+offset].setLinearCoef(v, col.coefsSRC[next_src]);
                        next_src++;
                    }
                }
                assert(next_rcc==RCCs.size());
                assert(next_src==SRCs.size());
                /*
                assert(col.coefsRCC.size()==RCCs.size());
                int offset=Config::FIX_VEHICLES*1+data->numNodes()-1;
                for (int i=0; i<RCCs.size(); ++i)
                    con[i+offset].setLinearCoef(v, col.coefsRCC[i]);
                */
                col.included=true;
            }
            solveBranchingFirst=false;
        }
        assert(columns.size()==colidxs.size());
        IloNumArray vals(env, colidxs.size());
        for (int i=0; i<colidxs.size(); ++i) {
            if (colidxs[i]==1)
                vals[i]=columns[i].expCost;
            else
                vals[i]=1e+6;
        }
        obj.setLinearCoefs(var, vals);
        vals.end();
        if (!cplex.solve()) {
            cerr<<"solveBranching(): failed to solve linear relax."<<endl;
            exit(-1);
        } else {
            Solution sol;
            sol.setValue(cplex.getObjValue());
            // saving optimal solution vector
            IloNumArray val(env);
            cplex.getValues(var, val);
            for (int i=0; i<val.getSize(); ++i)
                if (val[i]!=0)
                    sol.addRouteSolution(val[i], columns[i].route);
            val.end();
            return sol;
        }
    } catch (IloException& e) {
        cerr<<"solveBranching(): Concert exception caught: "<<e<<endl;
        exit(-1);
    } catch (...) {
        cerr<<"solveBranching(): unknown exception caught"<<endl;
        exit(-1);
    }
    cerr<<"solveBranching(): control should never reach here"<<endl;
    exit(-1);
    return Solution();        // to avoid return-type warnings
}

/* Solves the linear relaxation of the VRPSD by column generation. Computes the
 * exact linear relaxation. */
Solution LinearSolver::solveExact(bool sep3SRC, double upperBound,
        shared_ptr<IntegerSolver> isolver) {
    double lowerBound=0;
    const bool rootNode=arcsOn.size()==0&&arcsOff.size()==0;
    while (true) {
        currentSol=solveRMP();
        cout<<"e: current RMP sol value: "<<currentSol.value()<<endl;
        currentSol.printRoutes();
        if (Config::ENABLE_CUTS) {
            //cout<<"e: separating RCCs..."<<endl;
            SeparationRCC sep(currentSol);
            vector<RoundedCapacityCut> sepcuts=sep.solve();
            if (sepcuts.size()>0) {
                cout<<"e: "<<sepcuts.size()<<" RCCs found"<<endl;
                for (const auto& c : sepcuts)
                    addRCCut(c);
                continue;
            } else
                cout<<"e: no RCCs found"<<endl;
        }
        Pricing pricing(currentSol, arcsOn, arcsOff, RCCs, SRCs);
        vector<Route> routes=pricing.solveExact();
        if (routes.size()==0) {
            cout<<"e: no cols with neg rc found"<<endl;
            if (!sep3SRC || currentSol.value()>=upperBound)
                return currentSol;
            if (currentSol.value()>lowerBound&&rootNode) {
                lowerBound=currentSol.value();
                cout<<"LB: "<<lowerBound<<endl;
            }
            cout<<"e: separating 3SRCs..."<<endl;
            Separation3SRC sep(currentSol);
            vector<SubsetRowCut> sepcuts=sep.solve();
            if (sepcuts.size()>0) {
                cout<<"e: "<<sepcuts.size()<<" 3SRCs found"<<endl;
                vector<pair<SubsetRowCut,double>> srcs=findNBestCuts(sepcuts,8);
                double lbinc=(srcs[0].second-currentSol.value())
                        /currentSol.value();
                assert(lbinc>=0);
                cout<<"best SRC: "<<srcs[0].first<<"\t"<<(lbinc*100)
                        <<"\% increase in the LB"<<endl;
                if ((rootNode&&lbinc>Config::SRC_MIN_IMPROVEMENT)||
                        (!rootNode&&lbinc>4*Config::SRC_MIN_IMPROVEMENT)) {
                    cout<<"adding "<<srcs.size()<<" SRCs ..."<<endl;
                    for (const auto& src : srcs)
                        addSRCut(src.first);
                    continue;
                } else {
                    cout<<"very small LB increase by best SRC"<<endl;
                    return currentSol;
                }
            } else
                return currentSol;
        }
        cout<<"e: found "<<routes.size()<<" cols with neg rc:"<<endl;
        for (auto& r : routes) {
            cout<<r<<endl;
            addNewColumn(r);
        }
        if (isolver!=nullptr&&rootNode) {
            Solution isol=isolver->solve(columns);
            if (isol.feasible() && isol.value()<upperBound) {
                cout<<"new upper-bound obtained by the int solver"<<endl;
                cout<<"upper-bound improved: "<<isol.value()<<" (paths: "
                        <<isol.routes().size()<<")"<<endl;
                isol.print();
                upperBound=isol.value();
            }
            cout<<"gap (root node only!): "
                    <<100*((upperBound-lowerBound)/lowerBound)<<"\%"<<endl;
        }
    }
}

/* Solves the linear relaxation of the VRPSD by column generation.
 * This function solves the problem heuristically. Of course it does not
 * guarantee that the solution returned is optimal, nor that a feasible
 * solution is returned. If the returned solution is feasible then it is a
 * valid one (i.e., 'arcsOn' and 'arcsOff' are respected). If 'columns' has
 * size zero, then invalid routes are added (to achieve feasibility) and a
 * good set of initial routes is generated using the acceleration algorithm. */
Solution LinearSolver::solveHeuristic() {
    data->statIncPricingSolveHeuristic();
    if (columns.size()==0)
        genInitialColumns();
    int fresh=1;
    while (true) {
        currentSol=solveRMP();
        cout<<"h: current RMP sol value: "<<currentSol.value()<<endl;
        currentSol.printRoutes();
        if (Config::ENABLE_CUTS) {
            //cout<<"h: separating RCCs..."<<endl;
            SeparationRCC sep(currentSol);
            vector<RoundedCapacityCut> sepcuts=sep.solve();
            if (sepcuts.size()>0) {
                cout<<"h: "<<sepcuts.size()<<" RCCs found"<<endl;
                for (const auto& c : sepcuts)
                    addRCCut(c);
                continue;
            } else
                cout<<"h: no RCCs found"<<endl;
        }
        cout<<"h: solving pricing (h) from the RMP solution"<<endl;
        Pricing pricing(currentSol, arcsOn, arcsOff, RCCs, SRCs);
        vector<Route> routes=pricing.solveHeuristic();
        cout<<"h: found "<<routes.size()<<" cols with neg rc:"<<endl;
        for (auto& r : routes) {
            cout<<r<<endl;
            addNewColumn(r);
        }
        if (routes.size()<Config::COLS_LIMIT_HEUR) {
            fresh++;
            if (fresh==2)
                return currentSol=solveRMP();
        } else
            fresh=0;
    }
}

void LinearSolver::solvePRA() {
    // generate a good set of columns with the acceleration algorithm
    int maxLoad=data->maxLoad(true);
    int l=maxLoad*Config::CG_ACCEL_RATIO_MAXLOAD+0.5;
    int inc=max(1, (maxLoad-l)/Config::CG_ACCEL_ITERATIONS);
    while (l<maxLoad) {
        data->setACGMaxLoad(l);
        solveAccelerated();
        Solution sol=solveAccelerated();
        cout<<"accel. CG: done for l="<<l<<"\tsol. value="<<sol.value()<<endl;
        l+=inc;
    }
    data->setACGMaxLoad(data->maxLoad(true));
    cout<<"ACG finished"<<endl;
}

Solution LinearSolver::solveRMP(bool saveDuals) {
    cout<<"solveRMP(): RCCs: "<<RCCs.size()<<" SRCs: "<<SRCs.size()<<endl;
    //cout<<"solving the linear relaxation..."<<endl;
    try {
        /*
        for (auto& col : columns) {
            if (!col.included) {
                IloNumVar v(env, 0, IloInfinity);
                var.add(v);
                obj.setLinearCoef(v, col.compliesWithBranching?col.expCost
                        :1.1*col.expCost);
                if (Config::FIX_VEHICLES)
                    con[0].setLinearCoef(v, 1);
                for (int i=1; i<data->numNodes(); ++i)
                    if (col.route.visitsNode(i)>0)
                        con[Config::FIX_VEHICLES*1+i-1].setLinearCoef(
                                v, col.route.visitsNode(i));
                assert(col.coefsRCC.size()==RCCs.size());
                assert(col.coefsSRC.size()==SRCs.size());
                int offset=Config::FIX_VEHICLES*1+data->numNodes()-1;
                int next_rcc=0;
                int next_src=0;
                for (int i=0; i<RCCs.size()+SRCs.size(); ++i) {
                    assert(typeConstraint.at(i+offset)==CUT_RCC ||
                            typeConstraint[i+offset]==CUT_SRC);
                    if (typeConstraint[i+offset]==CUT_RCC) {
                        con[i+offset].setLinearCoef(v, col.coefsRCC[next_rcc]);
                        next_rcc++;
                    } else {
                        con[i+offset].setLinearCoef(v, col.coefsSRC[next_src]);
                        next_src++;
                    }
                }
                assert(next_rcc==RCCs.size());
                assert(next_src==SRCs.size());
                col.included=true;
            }
        }
        */
        // create variables and setup obj function coefficients
        IloNumVarArray objvars(env);
        IloNumArray objvals(env);
        for (const auto& col : columns) {
            if (!col.included) {
                IloNumVar v(env, 0, IloInfinity);
                var.add(v);
                objvars.add(v);
                objvals.add(col.compliesWithBranching?col.expCost
                        :1.1*col.expCost);
            }
        }
        if (objvars.getSize()>0) {
            obj.setLinearCoefs(objvars, objvals);
            if (Config::FIX_VEHICLES) {
                IloNumArray ones(env);
                for (int i=0; i<objvars.getSize(); ++i)
                    ones.add(1);
                con[0].setLinearCoefs(objvars, ones);
                ones.end();
            }
        } else {
            objvars.end();
            objvals.end();
            goto solve_directly;
        }
        // setup coefficients of the partitioning constraints
        typedef IloArray<IloNumVarArray> IloNumVarArray2;
        typedef IloArray<IloNumArray> IloNumArray2;
        {
            IloNumVarArray2 varspart(env);
            IloNumArray2 valspart(env);
            for (int i=0; i<data->numNodes(); ++i) {
                varspart.add(IloNumVarArray(env));
                valspart.add(IloNumArray(env));
            }
            for (int c=0; c<columns.size(); ++c) {
                const auto& col=columns[c];
                if (col.included)
                    continue;
                for (int i=1; i<data->numNodes(); ++i)
                    if (col.route.visitsNode(i)>0) {
                        varspart[i].add(var[c]);
                        valspart[i].add(col.route.visitsNode(i));
                    }
            }
            for (int i=0; i<data->numNodes(); ++i) {
                if (varspart[i].getSize()>0)
                    con[Config::FIX_VEHICLES*1+i-1].setLinearCoefs(
                            varspart[i], valspart[i]);
                varspart[i].end();
                valspart[i].end();
            }
            varspart.end();
            valspart.end();
        }
        // setup coefficients of cuts
        {
            IloNumVarArray2 varscuts(env);
            IloNumArray2 valscuts(env);
            for (int i=0; i<RCCs.size()+SRCs.size(); ++i) {
                varscuts.add(IloNumVarArray(env));
                valscuts.add(IloNumArray(env));
            }
            const int offset=Config::FIX_VEHICLES*1+data->numNodes()-1;
            for (int c=0; c<columns.size(); ++c) {
                auto& col=columns[c];
                if (col.included)
                    continue;
                assert(col.coefsRCC.size()==RCCs.size());
                assert(col.coefsSRC.size()==SRCs.size());
                col.included=true;
                int next_rcc=0;
                int next_src=0;
                for (int i=0; i<RCCs.size()+SRCs.size(); ++i) {
                    assert(typeConstraint.at(i+offset)==CUT_RCC ||
                            typeConstraint[i+offset]==CUT_SRC);
                    if (typeConstraint[i+offset]==CUT_RCC) {
                        varscuts[i].add(var[c]);
                        valscuts[i].add(col.coefsRCC[next_rcc]);
                        next_rcc++;
                    } else {
                        varscuts[i].add(var[c]);
                        valscuts[i].add(col.coefsSRC[next_src]);
                        next_src++;
                    }
                }
                assert(next_rcc==RCCs.size());
                assert(next_src==SRCs.size());
            }
            for (int i=0; i<RCCs.size()+SRCs.size(); ++i) {
                if (varscuts[i].getSize()>0)
                    con[i+offset].setLinearCoefs(varscuts[i], valscuts[i]);
                varscuts[i].end();
                valscuts[i].end();
            }
            varscuts.end();
            valscuts.end();
        }
        //cplex.exportModel("lp.lp");
solve_directly:
        if (!cplex.solve()) {
            cerr<<"solveRMP(): failed to solve linear relaxation"<<endl;
            exit(-1);
        } else {
            Solution sol;
            // saving optimal solution value
            sol.setValue(cplex.getObjValue());
            // saving optimal solution vector
            IloNumArray val(env);
            cplex.getValues(var, val);
            for (int i=0; i<val.getSize(); ++i)
                if (val[i]!=0) {
                    sol.addRouteSolution(val[i], columns[i].route);
                    columns[i].onceUsed=true;
                }
            assert(con.getSize()==Config::FIX_VEHICLES*1+data->numNodes()-1
                    +RCCs.size()+SRCs.size());
            // saving optimal solution dual variables
            if (saveDuals) {
                cplex.getDuals(val, con);
                if (Config::FIX_VEHICLES) {
                    for (int i=0; i<val.getSize(); ++i) {
                        if (i==0) {
                            assert(typeConstraint[0]==CTR_MAXV);
                            sol.addMaxVDual(val[i]);
                        } else if (i<data->numNodes()) {
                            assert(typeConstraint[i]==CTR_PTT);
                            sol.addPttDual(val[i]);
                        } else if (typeConstraint[i]==CUT_RCC)
                            sol.addRCCDual(val[i]);
                        else {
                            assert(typeConstraint[i]==CUT_SRC);
                            assert(val[i]<=0);
                            sol.addSRCDual(val[i]);
                        }
                    }
                } else {
                    for (int i=0; i<val.getSize(); ++i) {
                        if (i<data->numNodes()-1) {
                            assert(typeConstraint[i]==CTR_PTT);
                            sol.addPttDual(val[i]);
                        } else if (typeConstraint[i]==CUT_RCC)
                            sol.addRCCDual(val[i]);
                        else {
                            assert(typeConstraint[i]==CUT_SRC);
                            if (val[i]>1e-5)
                                cout<<"val["<<i<<"]: "<<val[i]<<endl;
                            assert(val[i]<=1e-5);
                            sol.addSRCDual(val[i]);
                        }
                    }
                }
            }
            val.end();
            // keep track of the best intermediate binary solution
            if (sol.feasible() && sol.binary() &&
                    (!intermBinary || sol.value()<intermBinarySol.value())) {
                cout<<"saving intermediate binary solution (solution value: "
                        <<sol.value()<<")"<<endl;
                intermBinarySol=sol;
                intermBinary=true;
            }
            return sol;
        }
    } catch (IloException& e) {
        cerr<<"solveRMP(): Concert exception caught: "<<e<<endl;
        exit(-1);
    } catch (...) {
        cerr<<"solveRMP(): unknown exception caught"<<endl;
        exit(-1);
    }
    cerr<<"solveRMP(): control should never reach here"<<endl;
    exit(-1);
    return Solution();        // to avoid return-type warnings
}

