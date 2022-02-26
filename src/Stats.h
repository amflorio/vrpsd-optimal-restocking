class Stats {
    private:
        bool acgFinished=false;
        int firstRule=0;
        int labelsExact=0;
        int labelsHeuristic=0;
        int pricingSolveBinary=0;
        int pricingSolveFeasible=0;
        int pricingSolveHeuristic=0;
        int pricingSolveLower=0;
        int zerothRuleExact=0;
        int zerothRuleHeuristic=0;
        int zerothRuleOPLRelax=0;
        int zerothRuleOPLRelaxTries=0;
    public:
        void add1stRule(int add) {firstRule+=add;}
        void inc0thRuleExact() {zerothRuleExact++;}
        void inc0thRuleHeuristic() {zerothRuleHeuristic++;}
        void inc0thRuleOPLRelax() {zerothRuleOPLRelax++;}
        void incLabelsExact() {labelsExact++;}
        void incLabelsHeuristic() {labelsHeuristic++;}
        void incOPLRelaxTries() {zerothRuleOPLRelaxTries++;}
        void incPricingSolveBinary() {pricingSolveBinary++;}
        void incPricingSolveFeasible() {pricingSolveFeasible++;}
        void incPricingSolveHeuristic() {pricingSolveHeuristic++;}
        void incPricingSolveLower() {pricingSolveLower++;}
        void print();
        void setACGFinished() {acgFinished=true;}
};

