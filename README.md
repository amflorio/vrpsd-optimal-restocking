# Elementary Branch-Price-and-Cut for the VRPSD under Optimal Restocking

## Description
This repository contains the source code and datasets to allow the replication of the results from the paper:

**Florio, A. M., Gendreau, M., Hartl, R. F., Minner, S., & Vidal, T. (2023). Recent Advances in Vehicle Routing with Stochastic Demands: Bayesian Learning for Correlated Demands and Elementary Branch-Price-and-Cut. European Journal of Operational Research 306(3), 1081-1093.**

The paper is also freely available at [https://arxiv.org/abs/2302.02538](https://arxiv.org/abs/2302.02538)

The code implements a branch-cut-and-price (BP&C) algorithm for the VRPSD under optimal restocking.

## Dependencies
The implementation requires:
* The [CPLEX Optimization Studio](https://www.ibm.com/ca-en/products/ilog-cplex-optimization-studio), which is free for academic use.

## Building and Running
A `Makefile` is provided for reference only. This should be adapted to match the specific host requirements, including `CPLEX` header files and libraries.

The code compiles into a single executable `vrpsd`. To facilitate scripting, a numerical code is assigned to each standard VRPSD instance. This mapping is available at the source file `Instance.cpp`. For example, instance `A-n32-k5` (Augerat, 1995) corresponds to the instance number `2001`.

To solve this instance, simply run:

```
$ ./vrpsd 0 2001
Instance info:
ID: 2001
N: 31
Q: 100
f: 1.000
round to integer: 0
max vehicles: 5
demand distribution: Poisson
expected demand: 410.000
load coeff.: 4.100
*** BRANCH AND PRICE ***
customer 1 demands in [4,39]
customer 2 demands in [5,42]
customer 3 demands in [0,19]
[...]
```

The app is output intensive, so we recommend redirecting the output to a file for posterior analysis. After a few minutes (even though instance `2001` has only 31 customers, it is a difficult instance as many SRCs need to be separated to close the optimality gap), the algorithm finishes and the optimal solution is printed alongside some statistics:

```
Solution value: 856.310
1.000 {0,12,1,16,30,0} n=4 l=72 ec=73.507 ap=73.486 rc=-0.273
1.000 {0,21,31,19,17,13,26,0} n=6 l=82 ec=157.816 ap=156.055 rc=-1.958
1.000 {0,20,5,25,10,15,29,27,0} n=7 l=91 ec=203.351 ap=195.048 rc=-4.344
1.000 {0,22,9,11,4,28,8,18,14,0} n=8 l=78 ec=238.666 ap=238.030 rc=-0.506
1.000 {0,7,2,3,23,6,24,0} n=6 l=87 ec=182.969 ap=179.151 rc=-0.335
sum of coefficients: 5.000
Accelerated Column Generation finished: 1
a-priori length of the DE solution: 787.808
Value of the Stochastic Solution (VSS) = 37.328	(4.177%)
BnP solution (value): 856.310
BnP solution (routes):
{0,12,1,16,30,0} n=4 l=72 ec=73.507 ap=73.486 rc=-0.273
{0,21,31,19,17,13,26,0} n=6 l=82 ec=157.816 ap=156.055 rc=-1.958
{0,20,5,25,10,15,29,27,0} n=7 l=91 ec=203.351 ap=195.048 rc=-4.344
{0,22,9,11,4,28,8,18,14,0} n=8 l=78 ec=238.666 ap=238.030 rc=-0.506
{0,7,2,3,23,6,24,0} n=6 l=87 ec=182.969 ap=179.151 rc=-0.335
verifying route expected costs with the SDP algorithm
Route: 73.507	SDP: 73.510	off by: -0.004%
(under Switch policy: 73.508	-0.003%)
Route: 157.816	SDP: 157.826	off by: -0.006%
(under Switch policy: 157.818	-0.005%)
Route: 203.351	SDP: 203.364	off by: -0.006%
(under Switch policy: 203.352	-0.006%)
Route: 238.666	SDP: 238.678	off by: -0.005%
(under Switch policy: 238.667	-0.004%)
Route: 182.969	SDP: 182.978	off by: -0.005%
(under Switch policy: 182.971	-0.004%)
Switch policy (total): -0.005%
TIME PROFILING
adding nodes: 151.003 secs
	... dcost: 132.054 secs
MEMORY PROFILING
profiling: peak memory allocated to labels: 0.782G
time elapsed: 577 secs
```

