#ifndef INTEGERSOLVER_H
#define INTEGERSOLVER_H

#include <vector>
#include "Data.h"
#include "LinearSolver.h"
#include "Solution.h"

ILOSTLBEGIN

struct Column;

class IntegerSolver {
    private:
        static std::shared_ptr<Data> data;
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloNumVarArray var;
        IloRangeArray con;
        IloObjective obj;
    public:
        IntegerSolver();
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        Solution solve(std::vector<Column>& cols);
};

#endif

