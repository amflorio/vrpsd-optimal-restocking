#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <ilcplex/ilocplex.h>
#include <memory>
#include <vector>
#include "Arc.h"
#include "IntegerSolver.h"
#include "RoundedCapacityCut.h"
#include "Route.h"
#include "Solution.h"
#include "SubsetRowCut.h"

ILOSTLBEGIN

struct Column {
    Route route;
    double expCost;
    bool included;
    bool isNew;
    bool inIntegerSolver=false;
    bool onceUsed=false;
    bool compliesWithBranching=true;
    std::vector<int> coefsRCC;
    std::vector<int> coefsSRC;
    Column(Route r, double c, bool i, bool n) : route{r}, expCost{c},
            included{i}, isNew{n} {}
    int addRCCCoef(const RoundedCapacityCut& cut) {
        int coef=0;
        for (const auto& a : cut.arcs)
            coef+=route.timesArc(a.i, a.j);
        coefsRCC.push_back(coef);
        return coef;
    }
    int addSRCCoef(const SubsetRowCut& cut) {
        int coef=0;
        assert(cut.nodes().size()==3);
        for (auto i : cut.nodes())
            coef+=route.visitsNode(i);
        coef/=2;
        coefsSRC.push_back(coef);
        return coef;
    }
};

class IntegerSolver;

class LinearSolver {
    private:
        const int CTR_MAXV=0;       // max vehicles constraint
        const int CTR_PTT=1;        // partitioning constraint
        const int CUT_RCC=2;        // rounded capacity cut
        const int CUT_SRC=3;        // subset row cut
        static std::shared_ptr<Data> data;
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloNumVarArray var;
        IloRangeArray con;
        std::vector<int> typeConstraint;
        IloObjective obj;
        std::vector<Column> columns;
        std::vector<RoundedCapacityCut>& RCCs;
        std::vector<SubsetRowCut>& SRCs;
        std::vector<Arc> arcsOn;
        std::vector<Arc> arcsOff;
        Solution currentSol;
        Solution intermBinarySol;
        bool intermBinary=false;
        bool solveBranchingFirst=true;
        void addSRCut(const SubsetRowCut& cut);
        void addRCCut(const RoundedCapacityCut& cut);
        void addNewColumn(Route route);
        std::vector<std::pair<SubsetRowCut, double>> findNBestCuts(
                const std::vector<SubsetRowCut>& cuts, int N) const;
        void genInitialColumns();
        Solution solveAccelerated();
    public:
        LinearSolver(std::vector<Column> cols,
                std::vector<RoundedCapacityCut>& vRCCs,
                std::vector<SubsetRowCut>& vSRCs, std::vector<Arc> aOn,
                std::vector<Arc> aOff);
        ~LinearSolver() {env.end();}
        Solution getIntermediateBinary() const {return intermBinarySol;}
        std::vector<Column> getNewColumns();
        bool hasIntermediateBinary() const {return intermBinary;}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        Solution solveBinary();
        Solution solveBranching(const std::vector<int>& colidxs);
        Solution solveExact(bool sep3SRC, double upperBound,
                std::shared_ptr<IntegerSolver> isolver);
        Solution solveHeuristic();
        void solvePRA();
        Solution solveRMP(bool saveDuals=true);
};

#endif

