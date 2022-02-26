#include "Config.h"

// Branch-and-bound settings
const double Config::INITIAL_UB=1024*1024;
const double Config::BRANCH_COEF_INT=0.05;
const double Config::EPS_BINARY_SOL=1e-3;

// Column generation settings
const double Config::INFEASIBLE_ROUTE_PENALTY=1024*1024;
const int Config::COLS_LIMIT_EXACT=20;
const int Config::COLS_LIMIT_HEUR=20;
const double Config::CG_ACCEL_RATIO_MAXLOAD=0.1;
const int Config::CG_ACCEL_ITERATIONS=16;
const int Config::THREADS_LSOLVER=1;

// Pricing settings
const double Config::MAX_MEM_LABELS=40;
const int Config::FSET_SIZE=4;

// Separation settings
const bool Config::ENABLE_CUTS=true;
const int Config::MAX_CUTS=2;
const double Config::EPS_EDGE_NONZERO=1e-4;
const bool Config::RCC_LEQ=true;

// Other settings
const double Config::MIN_TO_PRINT=1e-3;
const double Config::THRESH_PROB_ZERO=1e-5;     // set to 1e-5 to reproduce the
                                                // results from the paper
// Frequently changed settings
const bool Config::DETOUR_TO_DEPOT=false;
const bool Config::FIX_VEHICLES=true;
const double Config::LOAD_FACTOR=1.00;
const bool Config::HEURISTIC=false;
const double Config::SRC_MIN_IMPROVEMENT=0.000075;
const int Config::PROB_DIST=0;      // 0=Poisson;5=NegBinomial;6=MixedDists

