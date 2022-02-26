#include <exception>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include "Data.h"
#include "IntegerSolver.h"
#include "Solution.h"

using namespace std;

shared_ptr<Data> IntegerSolver::data;

ILOSTLBEGIN

IntegerSolver::IntegerSolver() : model{env}, cplex{env}, var{env}, con{env},
        obj{IloMinimize(env)} {
    //cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::Param::Parallel, -1);
    // rhs - max vehicles
    if (Config::FIX_VEHICLES)
        con.add(IloRange(env, -IloInfinity, data->maxVehicles()));
    // rhs - node visits
    for (int i=1; i<data->numNodes(); ++i)
        con.add(IloRange(env, 1, 1));           // set partitioning
    model.add(obj);
    model.add(con);
    cplex.extract(model);
}

Solution IntegerSolver::solve(vector<Column>& cols) {
    try {
        cout<<"IntegerSolver::solve(): adding new columns"<<endl;
        int added=0;
        for (auto& col : cols)
            if (!col.inIntegerSolver) {
                IloNumVar colVar(env, 0, 1, ILOINT);
                var.add(colVar);
                obj.setLinearCoef(colVar, col.expCost);
                if (Config::FIX_VEHICLES)
                    con[0].setLinearCoef(colVar, 1);
                for (int i=1; i<data->numNodes(); ++i)
                    con[Config::FIX_VEHICLES*1+i-1].setLinearCoef(
                            colVar, col.route.visitsNode(i));
                col.inIntegerSolver=true;
                added++;
            }
        cout<<"IntegerSolver::solve(): "<<added<<" columns added"<<endl;
        if (!cplex.solve()) {
            cerr<<"IntegerSolver::solve(): can't solve the model"<<endl;
            exit(-1);
        } else {
            Solution sol;
            sol.setValue(cplex.getObjValue());
            IloNumArray val(env);
            cplex.getValues(var, val);
            // ASSUMPTION: a column once enters 'cols' is never removed
            for (int i=0; i<val.getSize(); ++i)
                if (val[i]!=0)
                    sol.addRouteSolution(val[i], cols.at(i).route);
            return sol;
        }
    } catch (IloException& e) {
        cerr<<"Concert exception caught: "<<e<<endl;
        exit(-1);
    } catch (...) {
        cerr<<"Unknown exception caught"<<endl;
        exit(-1);
    }
    cerr<<"IntegerSolver::solve(): control should never reach here"<<endl;
    exit(-1);
    return Solution();        // to avoid return-type warnings
}

