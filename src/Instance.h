#ifndef INSTANCE_H
#define INSTANCE_H

#include <string>
#include <vector>

class Instance {
    private:
        const int PoissonDist=0;
        const int UniformDist=1;
        const int TwoPointDist=2;
        const int Louveaux3=3;
        const int Louveaux9=4;
        const int NegBinomialDist=5;
        const int MixedDists=6;
        const int XY=0;
        const int Explicit=1;
        const int FullMatrix=2;
        int TYPE=-1;
        bool roundInteger=false;
        double MY_ZERO=1e-12;
        double MY_ONE=1-MY_ZERO;
        double probZero;                        // for 2-point distribution
        int ID;
        int N;                                  // does NOT include depot
        int Q;
        double f;
        int maxV=0;
        std::vector<int> xcoords;
        std::vector<int> ycoords;
        std::vector<int> weights;
        std::vector<std::vector<double>> dmatrix;
        std::vector<int> dparams;
        std::vector<int> maxdemands;
        std::vector<std::vector<double>> demandprobs;
        std::vector<std::vector<int>> cvrpRoutes;
        int probDist;
        void checkData();
        void genDemandProbs();
        void genDistanceMatrix();
        void initialize();
        void loadCVRPOptimal(std::string filename);
        void loadEdgeFullMatrixInstance(std::string filename);
        void loadEdgeLowerRowInstance(std::string filename);
        void loadXYInstance(std::string filename);
    public:
        Instance(int i);
        void computeMaxVehicles();
        double dist(int i, int j) const {return dmatrix[i][j];}
        int expDemand(int i) const;
        int id() const {return ID;}
        std::vector<std::vector<int>> optCVRPRoutes() const {return cvrpRoutes;}
        int maxDemand(int i) const {return maxdemands[i];}
        int maxLoad() const {return (int)(Q*f+0.5);}
        int maxVehicles() const {return maxV;}
        int numCustomers() const {return N;}    // does NOT include depot
        void printDistanceMatrix() const;
        void printInfo() const;
        double prob(int i, int k) const;
        int vehicleQ() const {return Q;}
};

#endif

