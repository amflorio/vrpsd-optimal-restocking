#include <algorithm>
#include <chrono>
#include <random>
#include <vector>
#include "Config.h"
#include "Data.h"
#include "Knapsack.h"
#include "Label.h"
#include "Pricing.h"
#include "Route.h"

using namespace std;

shared_ptr<Data> Pricing::data;
PricingState Pricing::state;

Pricing::Pricing(const Solution& s, const vector<Arc>& aOn,
        const vector<Arc>& aOff, const vector<RoundedCapacityCut>& vRCCs,
        const vector<SubsetRowCut>& vSRCs) : sol{s}, arcsOn{aOn}, arcsOff{aOff},
        RCCs{vRCCs}, SRCs{vSRCs}, bnbForbidArc(data->numNodes(),
        vector<int>(data->numNodes(), 0)), arcDuals(data->numNodes(),
        vector<double>(data->numNodes(), 0)), maxArcDual(data->numNodes(), 0) {
    // setup bnbForbidArc according to the branching decisions
    for (const auto& a : arcsOn) {          // a.i must go to a.j, so ...
        for (int i=1; i<data->numNodes(); ++i)
            if (i!=a.i)
                bnbForbidArc[i][a.j]=1;     // noone else can go to a.j, and ...
        for (int j=1; j<data->numNodes(); ++j)
            if (j!=a.j)
                bnbForbidArc[a.i][j]=1;     // a.i can go nowhere else
    }
    for (const auto& a : arcsOff)
        bnbForbidArc[a.i][a.j]=1;
    // setup the coefficient of the arcs based on existing RCCs
    assert(RCCs.size()==sol.numRCCDuals());
    for (int i=0; i<RCCs.size(); ++i) {
        // sanity check
        assert((Config::RCC_LEQ && sol.getRCCDual(i)<1e-5)
                || (!Config::RCC_LEQ && sol.getRCCDual(i)>-1e-5));
        for (const auto& a : RCCs[i].arcs)
            arcDuals[a.i][a.j]+=sol.getRCCDual(i);
    }
    // when using >= RCCs, knapsack values have to be adjusted by the cut duals
    for (int i=0; i<data->numNodes(); ++i)
        for (int j=0; j<data->numNodes(); ++j)
            if (maxArcDual[i]<arcDuals[j][i])
                maxArcDual[i]=arcDuals[j][i];
    // RCSP bound
    genRCSPBound(Config::FSET_SIZE-(data->numNodes()>=100?2:0));
}

void Pricing::genRCSPBound(int fSize) {
    cout<<"generating strong RCSP bounds ... (F set size: "<<fSize<<")"<<endl;
    // first: select the customers to be in the F subset
    vector<pair<int, double>> ratios;
    for (int i=1; i<data->numNodes(); ++i)
        ratios.push_back(pair<int, double>(i,
                sol.getPttDual(i)/data->expDemand(i)));
    sort(ratios.begin(), ratios.end(),
            [] (const pair<int, double>& p1, const pair<int, double>& p2)
            {return p1.second>p2.second;});
    fset.clear();
    for (int i=0; i<fSize; ++i)
        fset.push_back(ratios[i].first);
    // second: compute and store the RCSP bounds for each subset of F
    rcsp.clear();
    for (int i=0; i<1<<fSize; ++i) {
        vector<int> forbidden(data->numNodes(), 0);
        for (int f=0; f<fset.size(); ++f)
            if (i>>f & 1)
                forbidden[fset[f]]=1;
        rcsp.push_back(rcspBounds(forbidden));
    }
}

// with 2-cycle elimination
vector<vector<double>> Pricing::rcspBounds(const vector<int>& forbidden) const {
    int L=data->maxLoad(true);
    double inf=numeric_limits<double>::infinity();
    vector<vector<double>> best(data->numNodes(), vector<double>(L+1, inf));
    vector<vector<int>> pred(data->numNodes(), vector<int>(L+1, 0));
    vector<vector<double>> secBest(data->numNodes(), vector<double>(L+1, inf));
    for (int q=0; q<=L; ++q)
        best[0][q]=0.0;
    for (int q=1; q<=L; ++q)
        for (int i=1; i<data->numNodes(); ++i) {
            best[i][q]=best[i][q-1];
            pred[i][q]=pred[i][q-1];
            secBest[i][q]=secBest[i][q-1];
            if (q<data->expDemand(i))
                continue;
            for (int j=0; j<data->numNodes(); ++j) {
                double v=inf;
                if (j!=i && !forbidden[j] && !bnbForbidArc[j][i]) {
                    if (pred[j][q-data->expDemand(i)]!=i) {     // no 2-cycle
                        v=best[j][q-data->expDemand(i)]+data->dist(j, i)
                                -arcDuals[j][i]-sol.getPttDual(i);
                    } else {        // best leads to cycle, pick 2nd best
                        v=secBest[j][q-data->expDemand(i)]+data->dist(j, i)
                                -arcDuals[j][i]-sol.getPttDual(i);
                    }
                }
                if (v<best[i][q]) {
                    // move best to secBest
                    secBest[i][q]=best[i][q];
                    // set new best
                    best[i][q]=v;
                    pred[i][q]=j;
                } else if (v<secBest[i][q]) {
                    // just need to update secBest
                    secBest[i][q]=v;
                }
            }
        }
    return best;
}

int Pricing::toRCSPIndex(const Label& l, int next) const {
    int idx=0;
    for (int i=0; i<fset.size(); ++i)
        if (l.visits(fset[i]) || fset[i]==next)
            idx+=1<<i;
    return idx;
}

double Pricing::knapsackBound(const Label& l) const {
    Knapsack kp;
    int rload=data->maxLoad(false)-l.load();
    kp.setCapacity(rload);
    for (int i=1; i<data->numNodes(); ++i)
        if (!l.visits(i) && data->expDemand(i)<=rload)
            kp.addItem(sol.getPttDual(i)+maxArcDual[i], data->expDemand(i));
    //return l.redcost()-kp.ub();
    return l.redcost()-kp.solve();
}

/* Attention: backward labelling. We cannot finish a route with a node that is
 * the starting node of an arc we branched on. This function is not efficiency
 * critical. */
void Pricing::addRootLabels(vector<vector<Label>>& labels) {
    vector<int> forbid(data->numNodes(), 0);
    for (auto& a : arcsOn)
        forbid[a.i]=1;
    for (int i=1; i<data->numNodes(); ++i) {
        if (!forbid[i]) {
            Label l(sol.getMaxVDual(), SRCs.size());
            l.addDepot();
            l.addNode(i, sol.getPttDual(i), arcDuals.at(i).at(0), SRCs);
            labels[i].push_back(move(l));
        }
    }
}

/* Returns true if and only if the partial route of Label l can be "closed",
 * i.e., finalized by adding the depot. We have to make sure there is no
 * pending arc on. We do not need to be obsessive with efficiency here, because
 * this is only called when the label has a negative reduced cost. */
bool Pricing::bnbAllowCloseExtension(const Label& l) const {
    // arcsOff do not play a role, should have been considered already
    for (auto& a : arcsOn)
        if (a.j==l.node())      // attention: backward labelling
            return false;       // cannot add 0 because have to add a.i
    return true;
}

double Pricing::estMemLabels(const vector<vector<Label>>& labels) const {
    long cap=0;
    long used=0;
    for (const auto& v : labels) {
        cap+=v.capacity();
        used+=v.size();
    }
    return (data->sizeUsedLabel()*used+data->sizeReservedLabel()*(cap-used))
            /(1024*1024*1024);
}

long Pricing::countLabels(const vector<vector<Label>>& labels) const {
    long cap=0;
    long used=0;
    for (const auto& v : labels) {
        cap+=v.capacity();
        used+=v.size();
    }
    data->memoryProfilingLabels(cap, used);
    return used;
}

vector<vector<Label>>& Pricing::initLabels(double newMaxVDual,
        const vector<double>& newPttDuals) {
    if (state.firstrun) {
        assert(state.labels.size()==0);
        state.firstrun=false;
        state.labels.insert(state.labels.end(), data->numNodes(),
                vector<Label>());
        addRootLabels(state.labels);
    } else {
        assert(state.labels.size()==data->numNodes());
        // recycling: update the cost of all labels saved in the state
        for (int i=1; i<data->numNodes(); ++i)
            for (auto& l : state.labels[i])
                l.updateCost(newMaxVDual, newPttDuals, arcDuals,
                        sol.getSRCDuals());
    }
    return state.labels;
}

void Pricing::printArcsOn() const {
    cout<<"arcs on:";
    for (auto& a : arcsOn)
        cout<<" "<<a;
    cout<<endl;
}

void Pricing::printArcsOff() const {
    cout<<"arcs off:";
    for (auto& a : arcsOff)
        cout<<" "<<a;
    cout<<endl;
}

vector<Route> Pricing::solveExact() {
    assert(data->maxLoad(true)==data->maxLoad(false));
    vector<Route> routes;
    cout<<"Pricing: solveExact(): initiating exhaustive search"<<endl;
    vector<vector<Label>> labels(data->numNodes(), vector<Label>());
    addRootLabels(labels);
    bool labelsLeft=true;
    vector<int> shuffled;
    for (int i=1; i<data->numNodes(); ++i)
        shuffled.push_back(i);
    std::random_shuffle(shuffled.begin(), shuffled.end());
    double rclim=0.0;
    double minrc=0.0;
    while (labelsLeft) {
        labelsLeft=false;
        for (int k=0; k<shuffled.size(); ++k) {
            if (k%3==0)
                countLabels(labels);
            int i=shuffled[k];
            if (routes.size()>0 && labels[i].size()>20000) {
                sort(labels[i].begin(), labels[i].end(),
                        [](const Label& l1, const Label& l2)
                        {return l1.redcost()<l2.redcost();});
                labels[i].erase(labels[i].begin()+10000, labels[i].end());
            }
            for (auto& l : labels[i]) {
                for (int kk=0; kk<shuffled.size(); ++kk) {
                    const int j=shuffled[kk];
                    if (l.visits(j) || bnbForbidArc[j][i] ||
                            l.load()+data->expDemand(j)>data->maxLoad(true))
                        continue;
                    int rload=data->maxLoad(true)-l.load();
                    const int Q=data->vehicleQ();
                    if (l.v(Q)+data->dist(j, i)-l.sumDuals()-arcDuals[j][i]
                            +rcspBound(toRCSPIndex(l, j), j, rload)>=rclim)
                        continue;
                    Label ext=l.extend(j, sol.getPttDual(j), arcDuals[j][i],
                            SRCs, sol.getSRCDuals());
                    if (ext.negativeCost() && bnbAllowCloseExtension(ext)) {
                        routes.push_back(ext.route());
                        if (routes.size()==Config::COLS_LIMIT_EXACT)
                            return routes;
                        if (routes.back().reducedCost()<minrc) {
                            minrc=routes.back().reducedCost();
                            rclim=min(rclim, minrc/4);
                        }
                    }
                    if (ext.v(Q)-ext.sumDuals()+sol.getPttDual(j)
                            +rcspBound(toRCSPIndex(ext, -1), j, rload)>=rclim)
                        continue;
                    if (knapsackBound(ext)>=rclim)
                        continue;
                    labels[j].push_back(move(ext));
                    labelsLeft=true;
                }
            }
            labels[i].clear();
        }
    }
    return routes;
}

vector<Route> Pricing::solveHeuristic() {
    vector<Route> routes;
    vector<vector<Label>>& labels=initLabels(sol.getMaxVDual(),
            sol.getPttDuals());
    bool labelsLeft=true;
    vector<int> shuffled;
    for (int i=1; i<data->numNodes(); ++i)
        shuffled.push_back(i);
    std::random_shuffle(shuffled.begin(), shuffled.end());
    double rclim=0.0;
    double minrc=0.0;
    while (labelsLeft) {
        labelsLeft=false;
        for (int k=0; k<shuffled.size(); ++k) {
            if (k%3==0)
                countLabels(labels);
            const int i=shuffled[k];
            if (labels[i].size()>20000) {
                sort(labels[i].begin(), labels[i].end(),
                        [](const Label& l1, const Label& l2)
                        {return l1.redcost()<l2.redcost();});
                labels[i].erase(labels[i].begin()+10000, labels[i].end());
            }
            for (auto& l : labels[i]) {
                if (l.extended())
                    continue;
                l.setExtended();
                for (int kk=0; kk<shuffled.size(); ++kk) {
                    int j=shuffled[kk];
                    if (l.visits(j) || bnbForbidArc[j][i] ||
                            l.load()+data->expDemand(j)>data->maxLoad(false))
                        continue;
                    int rload=data->maxLoad(false)-l.load();
                    const int Q=data->vehicleQ();
                    if (l.v(Q)+data->dist(j, i)-l.sumDuals()-arcDuals[j][i]
                            +rcspBound(toRCSPIndex(l, j), j, rload)>=rclim)
                        continue;
                    Label ext=l.extend(j, sol.getPttDual(j), arcDuals[j][i],
                            SRCs, sol.getSRCDuals());
                    if (ext.negativeCost() && bnbAllowCloseExtension(ext)) {
                        routes.push_back(ext.route());
                        if (routes.size()==Config::COLS_LIMIT_HEUR)
                            return routes;
                        if (routes.back().reducedCost()<minrc) {
                            minrc=routes.back().reducedCost();
                            rclim=minrc/4;
                        }
                        goto clear;
                    }
                    if (ext.v(Q)-ext.sumDuals()+sol.getPttDual(j)
                            +rcspBound(toRCSPIndex(ext, -1), j, rload)>=rclim)
                        continue;
                    if (knapsackBound(ext)>=rclim)
                        continue;
                    labels[j].push_back(move(ext));
                    labelsLeft=true;
                }
            }
clear:
            labels[i].clear();
        }
    }
    resetState();
    return routes;
}

