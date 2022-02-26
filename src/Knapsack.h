#include <vector>

struct KnapsackItem {
    double value;
    int weight;
    double density;
    KnapsackItem(double v, int w) : value{v}, weight{w}, density{v/w} {}
};

class Knapsack {
    private:
        std::vector<KnapsackItem> items;
        int C;              // capacity of the knapsack
    public:
        void addItem(double v, int w) {items.emplace_back(v, w);}
        void setCapacity(int cap) {C=cap;}
        double solve() const;
        double ub();
};

