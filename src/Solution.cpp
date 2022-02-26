#include <iostream>
#include "Config.h"
#include "Solution.h"

using namespace std;

shared_ptr<Data> Solution::data;

bool Solution::binary() const {
    if (val==Config::INITIAL_UB)
        return false;
    for (auto& rs : routesSol)
        if (rs.coef>=Config::EPS_BINARY_SOL
                && rs.coef<=1-Config::EPS_BINARY_SOL)
            return false;
    return true;
}

void Solution::print() const {
    printPttDuals();
    cout<<"Solution value: "<<val<<endl;
    printRoutes();
}

void Solution::printPttDuals() const {
    double sum=0.0;
    for (int i=1; i<dualsSP.size(); ++i) {
        double ratio=dualsSP[i]/data->expDemand(i);
        cout<<"lambda_"<<i<<" : "<<dualsSP[i]<<"\tratio: "<<ratio<<endl;
        sum+=dualsSP[i];
    }
    cout<<"sum: "<<sum<<endl;
    cout<<"maxVDual: "<<maxVDual<<endl;
}

void Solution::printRoutes() const {
    double coef=0.0;
    for (auto& rs : routesSol) {
        coef+=rs.coef;
        if (rs.route.expCost()==Config::INFEASIBLE_ROUTE_PENALTY) {
            cout<<"Warning: solution contains infeasible initial route"<<endl;
            cout<<rs.coef<<" "<<rs.route<<endl;
            cout<<"affecting solution value in: "
                    <<rs.coef*Config::INFEASIBLE_ROUTE_PENALTY<<endl;
        } else if (rs.coef>=Config::MIN_TO_PRINT)
            cout<<rs.coef<<" "<<rs.route<<endl;
    }
    cout<<"sum of coefficients: "<<coef<<endl;
}

vector<Route> Solution::routes() const {
    vector<Route> routes;
    for (auto& rs : routesSol)
        routes.push_back(rs.route);
    return routes;
}

bool Solution::feasible() {
    if (val==Config::INITIAL_UB)
        return false;
    for (auto& rs : routesSol)
        if (rs.route.expCost()==Config::INFEASIBLE_ROUTE_PENALTY)
            return false;
    return true;
}

