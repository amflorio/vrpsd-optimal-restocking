#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include "Config.h"
#include "Data.h"
#include "Label.h"

using namespace std;

shared_ptr<Data> Label::data;
const double Label::MY0=1e-6;              // got circularity with 1e-7 (?)

Label::Label(double maxVDual, int numSRCs) : _maxVDual{maxVDual},
        visitsSRC(numSRCs, 0), _visitsNode(data->numNodes(), 0) {
}

/* Base case: {node, 0} */
void Label::addNode(int node, double dual, double arcCoef,
        const vector<SubsetRowCut>& SRCs) {
    data->startProfilingAddNode();
    routeNodes.push_back(node);
    _visitsNode.at(node)++;
    _load=data->expDemand(node);
    _sumPttDuals=dual;
    _sumRCCDuals=arcCoef;
    apriori=data->dist(node, 0);
    vq=vector<double>(data->vehicleQ()+1, data->dist(node, 0));
    // update SRCs visits
    assert(visitsSRC.size()==SRCs.size());
    for (int k=0; k<SRCs.size(); ++k)
        if (SRCs[k].contains(node))
            visitsSRC[k]++;
    // update _redcost
    _redcost=data->dist(0, node);
    for (int k=data->minDemand(node); k<=data->maxDemand(node); ++k) {
        int rtrips=k<=data->vehicleQ()?0:
                (k<=2*data->vehicleQ()?1:
                (k<=3*data->vehicleQ()?2:3));
        _redcost+=data->prob(node, k)*(rtrips*(data->dist(node, 0)
                +data->dist(0, node))
                +vq.at(data->vehicleQ()*(rtrips+1)-k));
    }
    _redcost-=_sumPttDuals+_sumRCCDuals+_sumSRCDuals+_maxVDual;
    data->endProfilingAddNode();
}

/* General case: {node, next, ..., 0} */
void Label::addNode(int node, double dual, double arcCoef, const Label& next) {
    const int Q=data->vehicleQ();
    data->startProfilingAddNode();
    routeNodes.push_back(node);
    _visitsNode.at(node)++;
    _load+=data->expDemand(node);
    _sumPttDuals+=dual;
    _sumRCCDuals+=arcCoef;
    int nextnode=routeNodes.at(routeNodes.size()-2);
    apriori+=data->dist(node, nextnode);
    // update vq by performing one step of the SDP algorithm
    // indirect (replenishment at the depot) cost does not depend on q
    double vi=data->dist(node, 0)+data->dist(0, nextnode);
    if (Config::DETOUR_TO_DEPOT)
        vi=1024.0*1024.0;
    else {
        for (int k=data->minDemand(nextnode);k<=data->maxDemand(nextnode);++k) {
            int rtrips=k<=Q?0:(k<=2*Q?1:(k<=3*Q?2:3));
            vi+=data->prob(nextnode, k)*(rtrips*(data->dist(nextnode, 0)
                    +data->dist(0, nextnode))+next.v(Q*(rtrips+1)-k));
        }
    }
    data->startProfilingAddNodeDirectCost();
    // direct cost depends on q
    bool optimize=data->maxDemand(nextnode)<=Q;
    for (int q=Q; q>=0; --q) {
        double vd=data->dist(node, nextnode);
        if (optimize) {
            for (int k=data->minDemand(nextnode); k<=data->maxDemand(nextnode);
                    ++k) {
                int rtrips=k<=q?0:1;
                vd+=data->prob(nextnode, k)*(rtrips*(data->dist(nextnode, 0)
                        +data->dist(0, nextnode))+next.v(Q*rtrips+q-k));
            }
        } else {
            for (int k=data->minDemand(nextnode); k<=data->maxDemand(nextnode);
                    ++k) {
                int rtrips=k<=q?0:(k<=q+Q?1:(k<=q+2*Q?2:3));
                vd+=data->prob(nextnode, k)*(rtrips*(data->dist(nextnode, 0)
                        +data->dist(0, nextnode))+next.v(Q*rtrips+q-k));
            }
        }
        if (vd<vi)
            vq[q]=vd;
        else {                  // use threshold property
            while (q>=0) {
                vq[q]=vi;
                q--;
            }
            break;
        }
    }
    data->endProfilingAddNodeDirectCost();
    // update _redcost
    _redcost=data->dist(0, node);
    for (int k=data->minDemand(node); k<=data->maxDemand(node); ++k) {
        int rtrips=k<=Q?0:(k<=2*Q?1:(k<=3*Q?2:3));
        _redcost+=data->prob(node, k)*(rtrips*(data->dist(node, 0)
                +data->dist(0, node))+vq.at(Q*(rtrips+1)-k));
    }
    _redcost-=_sumPttDuals+_sumRCCDuals+_sumSRCDuals+_maxVDual;
    data->endProfilingAddNode();
}

Label Label::extend(int i, double dual_i, double arcCoef,
        const vector<SubsetRowCut>& SRCs, const vector<double>& SRCDuals)
        const {
    Label ext(*this);
    // update SRCs visits and duals - needs to be made BEFORE addNode
    for (int k=0; k<SRCs.size(); ++k)
        if (SRCs[k].contains(i)) {
            ext.visitsSRC[k]++;
            if (ext.visitsSRC[k]%2==0)
                ext._sumSRCDuals+=SRCDuals[k];
        }
    ext.addNode(i, dual_i, arcCoef, *this);
    ext.ext=false;  // important due to the position of setExtended in Pricing
    /*
    if (vq[data->vehicleQ()]+4e-3<apriori) {
        cout<<"vq["<<data->vehicleQ()<<"]: "<<vq[data->vehicleQ()]<<endl;
        cout<<"apriori: "<<apriori<<endl;
    }
    assert(vq[data->vehicleQ()]+4e-3>=apriori);
    */
    return ext;
}

Route Label::route() const {
    vector<int> routeRev=routeNodes;
    routeRev.push_back(0);
    std::reverse(std::begin(routeRev), std::end(routeRev));
    Route r(routeRev);
    r.setReducedCost(_redcost);
    r.setExpCost(_redcost+_sumPttDuals+_sumRCCDuals+_sumSRCDuals+_maxVDual);
    r.setLoad(_load);
    return r;
}

bool Label::spawns(const Route& r) const {
    vector<int> nodes=r.nodes();
    for (int i=0; i<routeNodes.size(); ++i)
        if (routeNodes.at(i)!=nodes.at(nodes.size()-1-i))
            return false;
    return true;
}

bool Label::spawns(const Route& r, int j) const {
    vector<int> nodes=r.nodes();
    int i=0;
    for (; i<routeNodes.size(); ++i)
        if (routeNodes.at(i)!=nodes.at(nodes.size()-1-i))
            return false;
    if (j!=nodes.at(nodes.size()-1-i))
        return false;
    return true;
}

void Label::updateCost(double newMaxVDual, const vector<double>& newPttDuals,
        const vector<vector<double>>& newCutArcCoeffs,
        const vector<double>& newSRCDuals) {
    assert(newSRCDuals.size()==visitsSRC.size());
    // TODO: can just keep the sum of the `old' duals
    double oldMaxVDual=_maxVDual;
    double oldPttDuals=_sumPttDuals;
    double oldCutDuals=_sumRCCDuals;
    double oldSRCDuals=_sumSRCDuals;
    _maxVDual=newMaxVDual;
    _sumPttDuals=0;
    for (int i=1; i<data->numNodes(); ++i)
        if (_visitsNode[i]>0)
            _sumPttDuals+=_visitsNode[i]*newPttDuals[i];
    _sumRCCDuals=0;
    // attention: routeNodes stores the route in reverse order
    for (int k=1; k<routeNodes.size(); ++k)
        _sumRCCDuals+=newCutArcCoeffs.at(routeNodes[k]).at(routeNodes[k-1]);
    _sumSRCDuals=0;
    for (int k=0; k<visitsSRC.size(); ++k)
        _sumSRCDuals+=(visitsSRC[k]/2)*newSRCDuals[k];
    _redcost=_redcost+oldPttDuals+oldCutDuals+oldSRCDuals+oldMaxVDual
            -_sumPttDuals-_sumRCCDuals-_sumSRCDuals-_maxVDual;
}

ostream& operator<<(ostream& os, const Label& l) {
    os<<"{["<<l.routeNodes.at(0);
    for (int i=1; i<l.routeNodes.size(); ++i)
        os<<"-"<<l.routeNodes.at(i);
    os<<"]}";
    return os;
}

