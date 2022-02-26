#ifndef ROUTE_H
#define ROUTE_H

#include <initializer_list>
#include <memory>
#include <vector>
#include "Arc.h"
#include "Data.h"

class Data;

/* Route
 * Stores information about a route that starts and finishes at the depot */
class Route {
    friend std::ostream& operator<<(std::ostream& os, const Route& r);
    private:
        static std::shared_ptr<Data> data;
        double _expCost=-1;
        int _load=-1;
        std::vector<int> _nodes;
        std::vector<int> vNode;
        double redCost=0;       // at the time the route is added to the RMP
        double apriori=0;
    public:
        Route();
        Route(std::initializer_list<char> route);
        Route(const std::vector<int>& route);
        void addNode(int node);
        double aPriori() const {return apriori;}
        bool complyWithArcOn(const Arc& a) const;
        bool complyWithArcOff(const Arc& a) const;
        double expCost() const;
        bool hasNode(int i) const {return vNode[i]>0;}
        int load() const {return _load;}
        std::vector<int> nodes() const {return _nodes;}
        int numVisits() const {return _nodes.size()-2;}
        double reducedCost() const {return redCost;}
        Route reverse() const;
        void setExpCost(double c) {_expCost=c;}
        void setLoad(int l) {_load=l;}
        void setReducedCost(double cost) {redCost=cost;}
        static void setData(const std::shared_ptr<Data> d);
        int timesArc(int i, int j) const;
        int visitsNode(int i) const {return vNode[i];}
        int operator[](int idx) const {return _nodes[idx];}
};

std::ostream& operator<<(std::ostream& os, const Route& r);

#endif

