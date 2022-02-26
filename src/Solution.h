#ifndef SOLUTION_H
#define SOLUTION_H

#include "Data.h"
#include "Route.h"

struct RouteSolution {
    double coef;
    Route route;
    RouteSolution(double c, Route r) : coef{c}, route{r} {}
};

class Solution {
    private:
        static std::shared_ptr<Data> data;
        std::vector<RouteSolution> routesSol;
        double val;
        double maxVDual=0.0;
        std::vector<double> dualsSP {-1};  // no constraint on # of vehicles
        std::vector<double> dualsRCC;
        std::vector<double> dualsSRC;
    public:
        void addMaxVDual(double dual) {maxVDual=dual;}
        void addPttDual(double dual) {dualsSP.push_back(dual);}
        void addRCCDual(double dual) {dualsRCC.push_back(dual);}
        void addSRCDual(double dual) {dualsSRC.push_back(dual);}
        void addRouteSolution(double coef, Route route)
            {routesSol.push_back(RouteSolution(coef, route));}
        bool binary() const;
        bool feasible();
        double getMaxVDual() const {return maxVDual;}
        double getRCCDual(int i) const {return dualsRCC.at(i);}
        double getPttDual(int i) const {return dualsSP.at(i);}
        const std::vector<double>& getPttDuals() const {return dualsSP;}
        const std::vector<double>& getSRCDuals() const {return dualsSRC;}
        int numRCCDuals() const {return dualsRCC.size();}
        void print() const;
        void printPttDuals() const;
        void printRoutes() const;
        std::vector<Route> routes() const;
        std::vector<RouteSolution> routesSolution() const {return routesSol;}
        static void setData(const std::shared_ptr<Data> d) {data=d;}
        void setValue(double v) {val=v;}
        double value() {return val;}
};

#endif

