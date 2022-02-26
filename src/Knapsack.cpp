#include <algorithm>
#include <iostream>
#include "Knapsack.h"

using namespace std;

double Knapsack::solve() const {
    vector<vector<double>> m(items.size()+1, vector<double>(C+1, 0.0));
    for (int i=1; i<=items.size(); ++i) {
        for (int q=0; q<=C; ++q) {
            int wi=items[i-1].weight;
            double vi=items[i-1].value;
            if (wi>q)
                m[i][q]=m[i-1][q];
            else
                m[i][q]=std::max(m[i-1][q], m[i-1][q-wi]+vi);
        }
    }
    return m[items.size()][C];
}

double Knapsack::ub() {
    sort(items.begin(), items.end(),
            [](const KnapsackItem& i1, const KnapsackItem& i2)
            {return i1.density>i2.density;});
    double ub=0;
    int qleft=C;        // remaining capacity
    for (int i=0; i<items.size(); ++i) {
        if (items[i].weight<=0) {
            cout<<"ub(): item weight must be positive"<<endl;
            exit(-1);
        }
        if (items[i].weight<=qleft) {
            ub+=items[i].value;
            qleft-=items[i].weight;
        } else {
            ub+=((1.0*qleft)/items[i].weight)*items[i].value;
            break;
        }
    }
    return ub;
}

