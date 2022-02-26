#ifndef CONFIG_H
#define CONFIG_H

class Config {
    public:
        // Heuristic toggle
        static const bool HEURISTIC;
        // Branch-and-bound settings
        static const double INITIAL_UB;
        static const double BRANCH_COEF_INT;
        static const double EPS_BINARY_SOL;
        // Column generation settings
        static const double INFEASIBLE_ROUTE_PENALTY;
        static const int COLS_LIMIT_EXACT;
        static const int COLS_LIMIT_HEUR;
        static const double CG_ACCEL_RATIO_MAXLOAD;
        static const int CG_ACCEL_ITERATIONS;
        static const int THREADS_LSOLVER;
        // Pricing settings
        static const double MAX_MEM_LABELS;
        static const int FSET_SIZE;
        // Separation settings
        static const bool ENABLE_CUTS;
        static const int MAX_CUTS;
        static const double EPS_EDGE_NONZERO;
        static const bool RCC_LEQ;
        static const double SRC_MIN_IMPROVEMENT;
        // Other settings
        static const bool DETOUR_TO_DEPOT;
        static const bool FIX_VEHICLES;
        static const double LOAD_FACTOR;
        static const double MIN_TO_PRINT;
        static const double THRESH_PROB_ZERO;
        static const int PROB_DIST;
};

#endif

