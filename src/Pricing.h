#ifndef PRICING_H
#define PRICING_H

#include <memory>
#include <vector>
#include "Arc.h"
#include "Data.h"
#include "Label.h"
#include "LinearSolver.h"
#include "Route.h"
#include "Solution.h"

struct PricingState {
    bool firstrun=true;
    std::vector<std::vector<Label>> labels;
};

class Pricing {
    private:
        static std::shared_ptr<Data> data;
        static PricingState state;
        Solution sol;
        std::vector<Arc> arcsOn;
        std::vector<Arc> arcsOff;
        std::vector<RoundedCapacityCut> RCCs;
        std::vector<SubsetRowCut> SRCs;
        std::vector<std::vector<int>> bnbForbidArc;
        std::vector<std::vector<std::vector<double>>> rcsp;
        std::vector<int> fset;
        std::vector<std::vector<double>> arcDuals;
        std::vector<double> maxArcDual;
        void addRootLabels(std::vector<std::vector<Label>>& labels);
        bool bnbAllowCloseExtension(const Label& l) const;
        bool bnbAllowExtension(const Label& l, int j) const;
        bool checkBound(int n, const Label& l) const;
        bool checkBound(int n, const Label& l, int j) const;
        long countLabels(const std::vector<std::vector<Label>>& labels) const;
        double estMemLabels(const std::vector<std::vector<Label>>& labels)
                const;
        void genRCSPBound(int fSize);
        std::vector<std::vector<Label>>& initLabels(double newMaxVDual,
                const std::vector<double>& newPttDuals);
        double knapsackBound(const Label& l) const;
        void printArcsOn() const;
        void printArcsOff() const;
        void processDominance(std::vector<std::vector<Label>>& labels,
                bool dualLoop, bool hdom) const;
        /*
        double rcspBound(int rcspIdx, int i, int l) const
                {return std::min(rcsp[rcspIdx][i][l], 0.0);}
        */
        double rcspBound(int rcspIdx, int i, int l) const
                {return rcsp[rcspIdx][i][l];}
        std::vector<std::vector<double>> rcspBounds(
                const std::vector<int>& forbidden) const;
        static void resetState() {state.firstrun=true; state.labels.clear();}
        int toRCSPIndex(const Label& l, int next) const;
    public:
        Pricing(const Solution& s, const std::vector<Arc>& aOn,
                const std::vector<Arc>& aOff,
                const std::vector<RoundedCapacityCut>& vRCCs,
                const std::vector<SubsetRowCut>& vSRCs);
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        std::vector<Route> solveExact();
        std::vector<Route> solveHeuristic();
};

#endif

