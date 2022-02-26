#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>
#include "Config.h"
#include "Data.h"
#include "Route.h"

using namespace std;

shared_ptr<Data> Route::data;

Route::Route() : vNode(data->numNodes(), 0) {
}

Route::Route(initializer_list<char> route) : vNode(data->numNodes(), 0) {
    for (int i : route)
        addNode(i);
}

Route::Route(const vector<int>& route) : vNode(data->numNodes(), 0) {
    for (int i=0; i<route.size(); ++i)
        addNode(route[i]);
}

void Route::addNode(int node) {
    if (_nodes.size()>0)
        apriori+=data->dist(_nodes.back(), node);
    _nodes.push_back(node);
    vNode[node]++;
    if (node!=0)
        _load+=data->expDemand(node);
}

bool Route::complyWithArcOn(const Arc& a) const {
    if (vNode[a.i]==0 && vNode[a.j]==0)
        return true;
    if ((vNode[a.i]==0 && vNode[a.j]>0) || (vNode[a.i]>0 && vNode[a.j]==0))
        return false;
    // both a.i and a.j are visited
    for (int n=1; n<_nodes.size(); ++n)
        if (_nodes[n-1]==a.i && _nodes[n]!=a.j)
            return false;
    return true;
}

bool Route::complyWithArcOff(const Arc& a) const {
    if (vNode[a.i]==0 || vNode[a.j]==0)
        return true;
    // both a.i and a.j are visited
    for (int n=1; n<_nodes.size(); ++n)
        if (_nodes[n-1]==a.i && _nodes[n]==a.j)
            return false;
    return true;
}

double Route::expCost() const {
    return _load>data->maxLoad(true)?Config::INFEASIBLE_ROUTE_PENALTY:_expCost;
}

Route Route::reverse() const {
    vector<int> routeRev=_nodes;
    std::reverse(std::begin(routeRev), std::end(routeRev));
    Route rev(routeRev);
    return rev;
}

void Route::setData(const shared_ptr<Data> d) {
    Route::data=d;
}

int Route::timesArc(int i, int j) const {
    if (vNode[i]==0 || vNode[j]==0)
        return 0;
    int times=0;
    for (int n=1; n<_nodes.size(); ++n)
        if (_nodes[n-1]==i && _nodes[n]==j)
            times++;
    return times;
}

ostream& operator<<(ostream& os, const Route& r) {
    double cost=r.expCost();
    if (cost==Config::INFEASIBLE_ROUTE_PENALTY)
        os<<"**";
    os<<"{"<<r._nodes[0];
    for (int i=1; i<r._nodes.size(); ++i)
        os<<","<<r._nodes[i];
    os<<"} n="<<r.numVisits()<<" l="<<r._load<<" ec="<<cost<<" ap="<<r.apriori
            <<" rc="<<r.redCost;
    return os;
}

