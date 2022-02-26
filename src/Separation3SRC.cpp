#include "Separation3SRC.h"

using namespace std;

shared_ptr<Data> Separation3SRC::data;

vector<SubsetRowCut> Separation3SRC::solve() const {
    vector<SubsetRowCut> cuts;
    for (int i=1; i<data->numNodes(); ++i)
        for (int j=i+1; j<data->numNodes(); ++j)
            for (int k=j+1; k<data->numNodes(); ++k)
                if (i!=j && i!=k && j!=k) {
                    double lhs=0.0;
                    for (const auto& rs : solution.routesSolution()) {
                        int sum_a=rs.route.visitsNode(i)+rs.route.visitsNode(j)
                                +rs.route.visitsNode(k);
                        lhs+=(sum_a/2)*rs.coef;
                    }
                    if (lhs>1+1e-3)
                        cuts.emplace_back(vector<int>({i,j,k}), 1, lhs-1);
                }
    /*
    sort(cuts.begin(), cuts.end(),
            [](const SubsetRowCut& c1, const SubsetRowCut& c2)
            {return c1.violation()>c2.violation();});
    */
    return cuts;
}

