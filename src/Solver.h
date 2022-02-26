#ifndef SOLVER_H
#define SOLVER_H

#include <memory>
#include "Arc.h"
#include "Data.h"
#include "LinearSolver.h"
#include "Solution.h"

struct BnBNode {
    Solution sol;
    std::vector<Arc> arcsOn;
    std::vector<Arc> arcsOff;
    BnBNode(const Solution& s, const std::vector<Arc>& aOn,
            const std::vector<Arc>& aOff) : sol{s}, arcsOn{aOn},
            arcsOff{aOff} {}
};

/* Solves the VRPSD with Branch-and-Bound */
class Solver {
    private:
        static std::shared_ptr<Data> data;
        std::vector<Column> columns;
        std::vector<RoundedCapacityCut> RCCs;
        std::vector<SubsetRowCut> SRCs;
        Solution upperBound;
        double lowerBound=0;
        std::vector<BnBNode> pending;
        void addNewColumnsCuts(const std::vector<Column>& newcols);
        bool assertCandidatesNotBranched(const std::vector<Arc>& candidates,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        void assertSolutionValid(const Solution& sol,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        void assertValidArcToBranch(const Arc& arc,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        void checkUpperBoundImproved(const LinearSolver& lsolver,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff);
        std::vector<Arc> arcToBranchCandidates(const Solution& sol);
        std::vector<Column> filterColsBranching(const std::vector<Column>& cols,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        std::vector<Column> filterCols(std::vector<Column> cols,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        std::vector<int> filterColsIndex(const std::vector<Column>& cols,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff) const;
        std::vector<Column> filterColsArcOn(const std::vector<Column>& cols,
                const Arc& arc) const;
        std::vector<Column> filterColsArcOff(const std::vector<Column>& cols,
                const Arc& arc) const;
        BnBNode nextNode();
        void printArcsOn(const std::vector<Arc>& arcsOn) const;
        void printArcsOff(const std::vector<Arc>& arcsOff) const;
        Arc selectBestCandidate(const std::vector<Arc>& candidates,
                const std::vector<Arc>& arcsOn,
                const std::vector<Arc>& arcsOff);
        void updateUpperBound(const Solution& newUB);
    public:
        Solver(const std::vector<Column>& cols, const Solution& ub)
                : columns{cols}, upperBound{ub} {}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        Solution solve();
};

#endif

