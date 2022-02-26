#ifndef LABEL_H
#define LABEL_H

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include "Data.h"
#include "Route.h"
#include "SubsetRowCut.h"

class Label {
    friend std::ostream& operator<<(std::ostream& os, const Label& l);
    friend inline bool operator<(const Label& l, const Label& r);
    private:
        static std::shared_ptr<Data> data;
        const static double MY0;
        int _load=0;
        double _redcost=0;              // assuming partial route is closed
        double _maxVDual=0.0;
        double _sumPttDuals=0;          // duals of the nodes visited
        double _sumRCCDuals=0;
        double _sumSRCDuals=0;
        std::vector<int> visitsSRC;
        double apriori=0;               // cost of the a priori tour
        std::vector<int> _visitsNode;
        std::vector<int> routeNodes;    // keeps track of the route
        std::vector<double> vq;
        bool ext=false;                 // whether label has been extended
        void addNode(int node, double dual, double arcCoef, const Label& next);
        int gamma(int k, int q) {return std::max(0,
                1+(int)std::floor(((k-q-1)*1.0)/data->vehicleQ()));}
    public:
        Label(double maxVDual, int numSRCs);
        void addDepot() {routeNodes.push_back(0);}
        void addNode(int node, double dual, double arcCoef,
                const std::vector<SubsetRowCut>& SRCs);
        double aPriori() const {return apriori;}
        Label extend(int i, double dual_i, double arcCoef,
                const std::vector<SubsetRowCut>& SRCs,
                const std::vector<double>& SRCDuals) const;
        bool extended() const {return ext;}
        int load() const {return _load;}
        bool negativeCost() const {return _redcost<0-MY0;}
        int node() const {return routeNodes.back();}
        double redcost() const {return _redcost;}
        Route route() const;
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        void setExtended() {ext=true;}
        bool spawns(const Route& r) const;
        bool spawns(const Route& r, int j) const;
        double sumDuals() const {return _sumPttDuals+_sumRCCDuals+_sumSRCDuals
                +_maxVDual;}
        void updateCost(double newMaxVDual, const std::vector<double>& newduals,
                const std::vector<std::vector<double>>& newCutArcCoeffs,
                const std::vector<double>& newSRCDuals);
        double v(int q) const {return vq[q];}
        bool visits(int j) const {return _visitsNode[j];}
};

std::ostream& operator<<(std::ostream& os, const Label& l);

#endif

