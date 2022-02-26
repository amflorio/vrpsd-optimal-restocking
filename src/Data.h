#ifndef DATA_H
#define DATA_H

#include <chrono>
#include <cmath>
#include <vector>
#include "Config.h"
#include "Instance.h"
#include "Route.h"
#include "Stats.h"

class Route;

class Data {
    private:
        Stats stats;         // store some statistics
        Instance instance;
        std::vector<int> maxdemands;
        std::vector<int> mindemands;
        std::vector<std::vector<double>> probs;
        std::vector<std::vector<double>> mindists;
        std::vector<std::vector<int>> mindistsplus;
        int maxLoadACG;
        int countBits(int s) const;
        void genDemands();
        void genMinDists();
        void genMinDistsPlus();
        // profiling
        std::chrono::steady_clock::time_point beginAddNode;
        std::chrono::steady_clock::time_point beginAddNodeDirectCost;
        long elapAddNode=0;
        long elapAddNodeDirectCost=0;
        double maxMemoryLabels=0.0;
        double sizeRLabel;
        double sizeULabel;
    public:
        Data(const Instance& i);
        int maxDemand(int i) const {return maxdemands[i];}
        int maxVehicles() const {return instance.maxVehicles();}
        int minDemand(int i) const {return mindemands[i];}
        int minVehicles(const std::vector<int>& nodes) const;
        int numNodes() const {return instance.numCustomers()+1;}
        double dist(int i, int j) const {return instance.dist(i, j);}
        int expDemand(int i) const {return instance.expDemand(i);}
        int getInstanceID() {return instance.id();}
        int maxLoad(bool strict) const {return strict?instance.maxLoad():
                std::min(instance.maxLoad(), maxLoadACG);}
        double minDistFrom(int i, int rload) const
                {return mindists.at(i).at(rload);}
        double minDistPlus(int i, int rload,
                const std::vector<int>& forbid) const;
        std::vector<Route> optCVRPRoutes() const;
        std::vector<Route> initialRoutes();
        void printDistanceMatrix() const {instance.printDistanceMatrix();}
        void printConfigParams() const;
        void printStats() {stats.print();}
        double prob(int i, int k) {return probs[i][k];}
        void setACGMaxLoad(int load) {maxLoadACG=load;}
        double sizeReservedLabel() const {return sizeRLabel;}
        double sizeUsedLabel() const {return sizeULabel;}
        void statIncPricingSolveBinary() {stats.incPricingSolveBinary();}
        void statIncPricingSolveFeasible() {stats.incPricingSolveFeasible();}
        void statIncPricingSolveHeuristic() {stats.incPricingSolveHeuristic();}
        void statIncPricingSolveLower() {stats.incPricingSolveLower();}
        void statSetACGFinished() {stats.setACGFinished();}
        int vehicleQ() const {return instance.vehicleQ();}
        // profiling functions
        void endProfilingAddNode();
        void endProfilingAddNodeDirectCost();
        void memoryProfilingLabels(long cap, long used);
        void printProfiling();
        void startProfilingAddNode();
        void startProfilingAddNodeDirectCost();
};

#endif

