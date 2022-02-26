#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "Instance.h"
#include "NegBinomial.h"
#include "Poisson.h"
#include "Route.h"

using namespace std;

Instance::Instance(int i) : ID{i} {
    f=Config::LOAD_FACTOR;
    probDist=Config::PROB_DIST; // default, unless otherwise specified
    TYPE=XY;                    // default
    if (i==1001)                // CMT
        loadXYInstance("../in/CMT1");
    else if (i==1002)
        loadXYInstance("../in/CMT2");
    else if (i==1003)
        loadXYInstance("../in/CMT3");
    else if (i==1004)
        loadXYInstance("../in/CMT4");
    else if (i==1005)
        loadXYInstance("../in/CMT5");
    else if (i==1011)
        loadXYInstance("../in/CMT11");
    else if (i==1012)
        loadXYInstance("../in/CMT12");
    else if (i==2001)           // Set A
        loadXYInstance("../in/A-n32-k5");
    else if (i==2002)
        loadXYInstance("../in/A-n33-k5");
    else if (i==2003)
        loadXYInstance("../in/A-n33-k6");
    else if (i==2004)
        loadXYInstance("../in/A-n34-k5");
    else if (i==2005)
        loadXYInstance("../in/A-n36-k5");
    else if (i==2006)
        loadXYInstance("../in/A-n37-k5");
    else if (i==2007)
        loadXYInstance("../in/A-n37-k6");
    else if (i==2008)
        loadXYInstance("../in/A-n38-k5");
    else if (i==2009)
        loadXYInstance("../in/A-n39-k5");
    else if (i==2010)
        loadXYInstance("../in/A-n39-k6");
    else if (i==2011)
        loadXYInstance("../in/A-n44-k6");
    else if (i==2012)
        loadXYInstance("../in/A-n45-k6");
    else if (i==2013)
        loadXYInstance("../in/A-n45-k7");
    else if (i==2014)
        loadXYInstance("../in/A-n46-k7");
    else if (i==2015)
        loadXYInstance("../in/A-n48-k7");
    else if (i==2016)
        loadXYInstance("../in/A-n53-k7");
    else if (i==2017)
        loadXYInstance("../in/A-n54-k7");
    else if (i==2018)
        loadXYInstance("../in/A-n55-k9");
    else if (i==2019)
        loadXYInstance("../in/A-n60-k9");
    else if (i==2020)
        loadXYInstance("../in/A-n61-k9");
    else if (i==2021)
        loadXYInstance("../in/A-n62-k8");
    else if (i==2022)
        loadXYInstance("../in/A-n63-k9");
    else if (i==2023)
        loadXYInstance("../in/A-n63-k10");
    else if (i==2024)
        loadXYInstance("../in/A-n64-k9");
    else if (i==2025)
        loadXYInstance("../in/A-n65-k9");
    else if (i==2026)
        loadXYInstance("../in/A-n69-k9");
    else if (i==2027)
        loadXYInstance("../in/A-n80-k10");
    else if (i==3001)           // Set B
        loadXYInstance("../in/B-n31-k5");
    else if (i==3002)
        loadXYInstance("../in/B-n34-k5");
    else if (i==3003)
        loadXYInstance("../in/B-n35-k5");
    else if (i==3004)
        loadXYInstance("../in/B-n38-k6");
    else if (i==3005)
        loadXYInstance("../in/B-n39-k5");
    else if (i==3006)
        loadXYInstance("../in/B-n41-k6");
    else if (i==3007)
        loadXYInstance("../in/B-n43-k6");
    else if (i==3008)
        loadXYInstance("../in/B-n44-k7");
    else if (i==3009)
        loadXYInstance("../in/B-n45-k5");
    else if (i==3010)
        loadXYInstance("../in/B-n45-k6");
    else if (i==3011)
        loadXYInstance("../in/B-n50-k7");
    else if (i==3012)
        loadXYInstance("../in/B-n50-k8");
    else if (i==3013)
        loadXYInstance("../in/B-n51-k7");
    else if (i==3014)
        loadXYInstance("../in/B-n52-k7");
    else if (i==3015)
        loadXYInstance("../in/B-n56-k7");
    else if (i==3016)
        loadXYInstance("../in/B-n57-k7");
    else if (i==3017)
        loadXYInstance("../in/B-n57-k9");
    else if (i==3018)
        loadXYInstance("../in/B-n63-k10");
    else if (i==3019)
        loadXYInstance("../in/B-n64-k9");
    else if (i==3020)
        loadXYInstance("../in/B-n66-k9");
    else if (i==3021)
        loadXYInstance("../in/B-n67-k10");
    else if (i==3022)
        loadXYInstance("../in/B-n68-k9");
    else if (i==3023)
        loadXYInstance("../in/B-n78-k10");
    else if (i==4001) {         // Louveaux instances
        loadXYInstance("../in/E031-09h-1");
        Q=84;
        probDist=Louveaux3;
    } else if (i==4002) {
        loadXYInstance("../in/E031-09h-2");
        Q=79;
        probDist=Louveaux3;
    } else if (i==4003) {
        loadXYInstance("../in/E031-09h-3");
        Q=84;
        probDist=Louveaux9;
    } else if (i==4004) {
        loadXYInstance("../in/E031-09h-4");
        Q=79;
        probDist=Louveaux9;
    } else if (i==4005) {
        loadXYInstance("../in/E031-09h-5");
        Q=59;
        probDist=Louveaux3;
    } else if (i==4006) {
        loadXYInstance("../in/E031-09h-6");
        Q=56;
        probDist=Louveaux3;
    } else if (i==4007) {
        loadXYInstance("../in/E031-09h-7");
        Q=59;
        probDist=Louveaux9;
    } else if (i==4008) {
        loadXYInstance("../in/E031-09h-8");
        Q=56;
        probDist=Louveaux9;
    } else if (i==4009) {
        loadXYInstance("../in/E051-05e-1");
        Q=139;
        probDist=Louveaux3;
    } else if (i==4010) {
        loadXYInstance("../in/E051-05e-2");
        Q=132;
        probDist=Louveaux3;
    } else if (i==4011) {
        loadXYInstance("../in/E051-05e-3");
        Q=139;
        probDist=Louveaux9;
    } else if (i==4012) {
        loadXYInstance("../in/E051-05e-4");
        Q=132;
        probDist=Louveaux9;
    } else if (i==4013) {
        loadXYInstance("../in/E051-05e-5");
        Q=99;
        probDist=Louveaux3;
    } else if (i==4014) {
        loadXYInstance("../in/E051-05e-6");
        Q=93;
        probDist=Louveaux3;
    } else if (i==4015) {
        loadXYInstance("../in/E051-05e-7");
        Q=99;
        probDist=Louveaux9;
    } else if (i==4016) {
        loadXYInstance("../in/E051-05e-8");
        Q=93;
        probDist=Louveaux9;
    } else if (i==4017) {
        loadXYInstance("../in/E076-07s");
        Q=209;
        probDist=Louveaux3;
    } else if (i==4018) {
        loadXYInstance("../in/E076-07s");
        Q=198;
        probDist=Louveaux3;
    } else if (i==4019) {
        loadXYInstance("../in/E076-07s");
        Q=209;
        probDist=Louveaux9;
    } else if (i==4020) {
        loadXYInstance("../in/E076-07s");
        Q=198;
        probDist=Louveaux9;
    } else if (i==4021) {
        loadXYInstance("../in/E076-07s");
        Q=148;
        probDist=Louveaux3;
    } else if (i==4022) {
        loadXYInstance("../in/E076-07s");
        Q=139;
        probDist=Louveaux3;
    } else if (i==4023) {
        loadXYInstance("../in/E076-07s");
        Q=148;
        probDist=Louveaux9;
    } else if (i==4024) {
        loadXYInstance("../in/E076-07s");
        Q=139;
        probDist=Louveaux9;
    } else if (i==4025) {
        loadXYInstance("../in/E101-08e");
        Q=278;
        probDist=Louveaux3;
    } else if (i==4026) {
        loadXYInstance("../in/E101-08e");
        Q=264;
        probDist=Louveaux3;
    } else if (i==4027) {
        loadXYInstance("../in/E101-08e");
        Q=278;
        probDist=Louveaux9;
    } else if (i==4028) {
        loadXYInstance("../in/E101-08e");
        Q=264;
        probDist=Louveaux9;
    } else if (i==4029) {
        loadXYInstance("../in/E101-08e");
        Q=197;
        probDist=Louveaux3;
    } else if (i==4030) {
        loadXYInstance("../in/E101-08e");
        Q=186;
        probDist=Louveaux3;
    } else if (i==4031) {
        loadXYInstance("../in/E101-08e");
        Q=197;
        probDist=Louveaux9;
    } else if (i==4032) {
        loadXYInstance("../in/E101-08e");
        Q=186;
        probDist=Louveaux9;
    } else if (i==4033) {
        loadEdgeFullMatrixInstance("../in/A034-02f-1");
        Q=92;
        probDist=Louveaux3;
    } else if (i==4034) {
        loadEdgeFullMatrixInstance("../in/A034-02f-2");
        Q=87;
        probDist=Louveaux3;
    } else if (i==4035) {
        loadEdgeFullMatrixInstance("../in/A034-02f-3");
        Q=92;
        probDist=Louveaux9;
    } else if (i==4036) {
        loadEdgeFullMatrixInstance("../in/A034-02f-4");
        Q=87;
        probDist=Louveaux9;
    } else if (i==4037) {
        loadEdgeFullMatrixInstance("../in/A034-02f-5");
        Q=65;
        probDist=Louveaux3;
    } else if (i==4038) {
        loadEdgeFullMatrixInstance("../in/A034-02f-6");
        Q=62;
        probDist=Louveaux3;
    } else if (i==4039) {
        loadEdgeFullMatrixInstance("../in/A034-02f-7");
        Q=65;
        probDist=Louveaux9;
    } else if (i==4040) {
        loadEdgeFullMatrixInstance("../in/A034-02f-8");
        Q=62;
        probDist=Louveaux9;
    } else if (i==4041) {
        loadEdgeFullMatrixInstance("../in/A048-03f-1");
        Q=131;
        probDist=Louveaux3;
    } else if (i==4042) {
        loadEdgeFullMatrixInstance("../in/A048-03f-2");
        Q=124;
        probDist=Louveaux3;
    } else if (i==4043) {
        loadEdgeFullMatrixInstance("../in/A048-03f-3");
        Q=131;
        probDist=Louveaux9;
    } else if (i==4044) {
        loadEdgeFullMatrixInstance("../in/A048-03f-4");
        Q=124;
        probDist=Louveaux9;
    } else if (i==4045) {
        loadEdgeFullMatrixInstance("../in/A048-03f-5");
        Q=93;
        probDist=Louveaux3;
    } else if (i==4046) {
        loadEdgeFullMatrixInstance("../in/A048-03f-6");
        Q=88;
        probDist=Louveaux3;
    } else if (i==4047) {
        loadEdgeFullMatrixInstance("../in/A048-03f-7");
        Q=93;
        probDist=Louveaux9;
    } else if (i==4048) {
        loadEdgeFullMatrixInstance("../in/A048-03f-8");
        Q=88;
        probDist=Louveaux9;
    } else if (i==5001)         // Set E
        loadEdgeLowerRowInstance("../in/E-n13-k4");
    else if (i==5002)
        loadXYInstance("../in/E-n22-k4");
    else if (i==5003) {
        // we do not experiment with this instance because of numerical issues
        // last run: expected cost was much smaller than apriori cost
        cerr<<"instance 5003: fix numerical issues before running"<<endl;
        exit(-1);
        MY_ONE/=1e+1;
        loadXYInstance("../in/E-n23-k3");
    } else if (i==5004) {
        // we do not experiment with this instance because of numerical issues
        // last run: expected cost was much smaller than apriori cost
        cerr<<"instance 5004: fix numerical issues before running"<<endl;
        exit(-1);
        MY_ONE/=1e+1;
        loadXYInstance("../in/E-n30-k3");
    } else if (i==5005)
        loadEdgeLowerRowInstance("../in/E-n31-k7");
    else if (i==5006)
        loadXYInstance("../in/E-n33-k4");
    else if (i==5007)
        loadXYInstance("../in/E-n51-k5");
    else if (i==5008)
        loadXYInstance("../in/E-n76-k7");
    else if (i==5009)
        loadXYInstance("../in/E-n76-k8");
    else if (i==5010)
        loadXYInstance("../in/E-n76-k10");
    else if (i==5011)
        loadXYInstance("../in/E-n76-k14");
    else if (i==5012)
        loadXYInstance("../in/E-n101-k8");
    else if (i==5013)
        loadXYInstance("../in/E-n101-k14");
    else if (i==6001)           // Set P
        loadXYInstance("../in/P-n16-k8");
    else if (i==6002)
        loadXYInstance("../in/P-n19-k2");
    else if (i==6003)
        loadXYInstance("../in/P-n20-k2");
    else if (i==6004)
        loadXYInstance("../in/P-n21-k2");
    else if (i==6005)
        loadXYInstance("../in/P-n22-k2");
    else if (i==6006)
        loadXYInstance("../in/P-n22-k8");
    else if (i==6007)
        loadXYInstance("../in/P-n23-k8");
    else if (i==6008)
        loadXYInstance("../in/P-n40-k5");
    else if (i==6009)
        loadXYInstance("../in/P-n45-k5");
    else if (i==6010)
        loadXYInstance("../in/P-n50-k7");
    else if (i==6011)
        loadXYInstance("../in/P-n50-k8");
    else if (i==6012)
        loadXYInstance("../in/P-n50-k10");
    else if (i==6013)
        loadXYInstance("../in/P-n51-k10");
    else if (i==6014)
        loadXYInstance("../in/P-n55-k7");
    else if (i==6015)
        loadXYInstance("../in/P-n55-k8");
    else if (i==6016)
        loadXYInstance("../in/P-n55-k10");
    else if (i==6017)
        loadXYInstance("../in/P-n55-k15");
    else if (i==6018)
        loadXYInstance("../in/P-n60-k10");
    else if (i==6019)
        loadXYInstance("../in/P-n60-k15");
    else if (i==6020)
        loadXYInstance("../in/P-n65-k10");
    else if (i==6021)
        loadXYInstance("../in/P-n70-k10");
    else if (i==6022)
        loadXYInstance("../in/P-n76-k4");
    else if (i==6023)
        loadXYInstance("../in/P-n76-k5");
    else if (i==6024)
        loadXYInstance("../in/P-n101-k4");
    else if (i==7001)           // Set A with Q set to 80
        loadXYInstance("../in/A-n32-k5-D5050");
    else if (i==7002)
        loadXYInstance("../in/A-n33-k5-D5050");
    else if (i==7003)
        loadXYInstance("../in/A-n33-k6-D5050");
    else if (i==7004)
        loadXYInstance("../in/A-n34-k5-D5050");
    else if (i==7005)
        loadXYInstance("../in/A-n36-k5-D5050");
    else if (i==7006)
        loadXYInstance("../in/A-n37-k5-D5050");
    else if (i==7007)
        loadXYInstance("../in/A-n37-k6-D5050");
    else if (i==7008)
        loadXYInstance("../in/A-n38-k5-D5050");
    else if (i==7009)
        loadXYInstance("../in/A-n39-k5-D5050");
    else if (i==7010)
        loadXYInstance("../in/A-n39-k6-D5050");
    else if (i==7011)
        loadXYInstance("../in/A-n44-k6-D5050");
    else if (i==7012)
        loadXYInstance("../in/A-n45-k6-D5050");
    else if (i==7013)
        loadXYInstance("../in/A-n45-k7-D5050");
    else if (i==7014)
        loadXYInstance("../in/A-n46-k7-D5050");
    else if (i==7015)
        loadXYInstance("../in/A-n48-k7-D5050");
    else if (i==7016)
        loadXYInstance("../in/A-n53-k7-D5050");
    else if (i==7017)
        loadXYInstance("../in/A-n54-k7-D5050");
    else if (i==7018)
        loadXYInstance("../in/A-n55-k9-D5050");
    else if (i==7019)
        loadXYInstance("../in/A-n60-k9-D5050");
    else if (i==7020)
        loadXYInstance("../in/A-n61-k9-D5050");
    else if (i==7021)
        loadXYInstance("../in/A-n62-k8-D5050");
    else if (i==7022)
        loadXYInstance("../in/A-n63-k9-D5050");
    else if (i==7023)
        loadXYInstance("../in/A-n63-k10-D5050");
    else if (i==7024)
        loadXYInstance("../in/A-n64-k9-D5050");
    else if (i==7025)
        loadXYInstance("../in/A-n65-k9-D5050");
    else if (i==7026)
        loadXYInstance("../in/A-n69-k9-D5050");
    else if (i==7027)
        loadXYInstance("../in/A-n80-k10-D5050");
    else if (i==8001)           // First Solomon instance
        loadXYInstance("../in/R101.25");
    else if (i==8002)           // 8001 with Q set to 180
        loadXYInstance("../in/R101.25-Q180");
    else if (i==8003)           // Q set to 160
        loadXYInstance("../in/R101.25-Q160");
    else if (i==8004)           // Q set to 140
        loadXYInstance("../in/R101.25-Q140");
    else if (i==8005)           // Q set to 120
        loadXYInstance("../in/R101.25-Q120");
    else if (i==8006)           // Q set to 100
        loadXYInstance("../in/R101.25-Q100");
    else if (i==8007)           // Q set to 80
        loadXYInstance("../in/R101.25-Q080");
    else if (i==8008)           // Q set to 60
        loadXYInstance("../in/R101.25-Q060");
    else if (i==8009)           // Q set to 40
        loadXYInstance("../in/R101.25-Q040");
    else if (i==9001) {         // Louveaux instances
        loadXYInstance("../in/E031-09h-1");
        Q=84;
    } else if (i==9002) {
        loadXYInstance("../in/E031-09h-2");
        Q=79;
    } else if (i==9003) {
        loadXYInstance("../in/E031-09h-3");
        Q=84;
    } else if (i==9004) {
        loadXYInstance("../in/E031-09h-4");
        Q=79;
    } else if (i==9005) {
        loadXYInstance("../in/E031-09h-5");
        Q=59;
    } else if (i==9006) {
        loadXYInstance("../in/E031-09h-6");
        Q=56;
    } else if (i==9007) {
        loadXYInstance("../in/E031-09h-7");
        Q=59;
    } else if (i==9008) {
        loadXYInstance("../in/E031-09h-8");
        Q=56;
    } else if (i==9009) {
        loadXYInstance("../in/E051-05e-1");
        Q=139;
    } else if (i==9010) {
        loadXYInstance("../in/E051-05e-2");
        Q=132;
    } else if (i==9011) {
        loadXYInstance("../in/E051-05e-3");
        Q=139;
    } else if (i==9012) {
        loadXYInstance("../in/E051-05e-4");
        Q=132;
    } else if (i==9013) {
        loadXYInstance("../in/E051-05e-5");
        Q=99;
    } else if (i==9014) {
        loadXYInstance("../in/E051-05e-6");
        Q=93;
    } else if (i==9015) {
        loadXYInstance("../in/E051-05e-7");
        Q=99;
    } else if (i==9016) {
        loadXYInstance("../in/E051-05e-8");
        Q=93;
    } else if (i==9017) {
        loadXYInstance("../in/E076-07s");
        Q=209;
    } else if (i==9018) {
        loadXYInstance("../in/E076-07s");
        Q=198;
    } else if (i==9019) {
        loadXYInstance("../in/E076-07s");
        Q=209;
    } else if (i==9020) {
        loadXYInstance("../in/E076-07s");
        Q=198;
    } else if (i==9021) {
        loadXYInstance("../in/E076-07s");
        Q=148;
    } else if (i==9022) {
        loadXYInstance("../in/E076-07s");
        Q=139;
    } else if (i==9023) {
        loadXYInstance("../in/E076-07s");
        Q=148;
    } else if (i==9024) {
        loadXYInstance("../in/E076-07s");
        Q=139;
    } else if (i==9025) {
        loadXYInstance("../in/E101-08e");
        Q=278;
    } else if (i==9026) {
        loadXYInstance("../in/E101-08e");
        Q=264;
    } else if (i==9027) {
        loadXYInstance("../in/E101-08e");
        Q=278;
    } else if (i==9028) {
        loadXYInstance("../in/E101-08e");
        Q=264;
    } else if (i==9029) {
        loadXYInstance("../in/E101-08e");
        Q=197;
    } else if (i==9030) {
        loadXYInstance("../in/E101-08e");
        Q=186;
    } else if (i==9031) {
        loadXYInstance("../in/E101-08e");
        Q=197;
    } else if (i==9032) {
        loadXYInstance("../in/E101-08e");
        Q=186;
    } else if (i==9033) {
        loadEdgeFullMatrixInstance("../in/A034-02f-1");
        Q=92;
    } else if (i==9034) {
        loadEdgeFullMatrixInstance("../in/A034-02f-2");
        Q=87;
    } else if (i==9035) {
        loadEdgeFullMatrixInstance("../in/A034-02f-3");
        Q=92;
    } else if (i==9036) {
        loadEdgeFullMatrixInstance("../in/A034-02f-4");
        Q=87;
    } else if (i==9037) {
        loadEdgeFullMatrixInstance("../in/A034-02f-5");
        Q=65;
    } else if (i==9038) {
        loadEdgeFullMatrixInstance("../in/A034-02f-6");
        Q=62;
    } else if (i==9039) {
        loadEdgeFullMatrixInstance("../in/A034-02f-7");
        Q=65;
    } else if (i==9040) {
        loadEdgeFullMatrixInstance("../in/A034-02f-8");
        Q=62;
    } else if (i==9041) {
        loadEdgeFullMatrixInstance("../in/A048-03f-1");
        Q=131;
    } else if (i==9042) {
        loadEdgeFullMatrixInstance("../in/A048-03f-2");
        Q=124;
    } else if (i==9043) {
        loadEdgeFullMatrixInstance("../in/A048-03f-3");
        Q=131;
    } else if (i==9044) {
        loadEdgeFullMatrixInstance("../in/A048-03f-4");
        Q=124;
    } else if (i==9045) {
        loadEdgeFullMatrixInstance("../in/A048-03f-5");
        Q=93;
    } else if (i==9046) {
        loadEdgeFullMatrixInstance("../in/A048-03f-6");
        Q=88;
    } else if (i==9047) {
        loadEdgeFullMatrixInstance("../in/A048-03f-7");
        Q=93;
    } else if (i==9048) {
        loadEdgeFullMatrixInstance("../in/A048-03f-8");
        Q=88;
    } else if (i==10001)          // Set X
        loadXYInstance("../in/X-n101-k25");
    else if (i==10002)
        loadXYInstance("../in/X-n106-k14");
    else if (i==10003)
        loadXYInstance("../in/X-n110-k13");
    else if (i==10004)
        loadXYInstance("../in/X-n115-k10");
    else if (i==10005)
        loadXYInstance("../in/X-n120-k6");
    else if (i==10006)
        loadXYInstance("../in/X-n125-k30");
    else if (i==10007)
        loadXYInstance("../in/X-n129-k18");
    else if (i==10008)
        loadXYInstance("../in/X-n134-k13");
    else if (i==10009)
        loadXYInstance("../in/X-n139-k10");
    else if (i==10010)
        loadXYInstance("../in/X-n143-k7");
    else if (i==10011)
        loadXYInstance("../in/X-n148-k46");
    else if (i==10012)
        loadXYInstance("../in/X-n153-k22");
    else if (i==10013)
        loadXYInstance("../in/X-n157-k13");
    else if (i==10014)
        loadXYInstance("../in/X-n162-k11");
    else if (i==10015)
        loadXYInstance("../in/X-n167-k10");
    else if (i==10016)
        loadXYInstance("../in/X-n172-k51");
    else if (i==10017)
        loadXYInstance("../in/X-n176-k26");
    else if (i==10018)
        loadXYInstance("../in/X-n181-k23");
    else if (i==10019)
        loadXYInstance("../in/X-n186-k15");
    else if (i==10020)
        loadXYInstance("../in/X-n190-k8");
    else if (i==10021)
        loadXYInstance("../in/X-n195-k51");
    else if (i==10022)
        loadXYInstance("../in/X-n200-k36");
    else if (i==10023)
        loadXYInstance("../in/X-n204-k19");
    else if (i==10024)
        loadXYInstance("../in/X-n209-k16");
    else if (i==10025)
        loadXYInstance("../in/X-n214-k11");
    else if (i==10026)
        loadXYInstance("../in/X-n219-k73");
    else if (i==10027)
        loadXYInstance("../in/X-n223-k34");
    else if (i==10031)
        loadXYInstance("../in/X-n242-k48");
    else {
        cerr<<"instance "<<i<<" not implemented"<<endl;
        exit(-1);
    }
    Q=(Q/Config::LOAD_FACTOR)+0.5;
    //Q=(0.75*(Q/Config::LOAD_FACTOR))+0.5;
    if ((i>=4001 && i<=4999) || (i>=9001 && i<=9999))
        roundInteger=true;
    initialize();
}

void Instance::computeMaxVehicles() {
    double sumD=0.0;
    for (int i=1; i<=N; ++i)
        sumD+=expDemand(i);
    maxV=ceil(sumD/(f*Q));
}

void Instance::checkData() {
    if (TYPE==XY) {
        if (xcoords.size()!=N+1||ycoords.size()!=N+1) {
            cerr<<"data inconsistency: xcoords or ycoords wrong size"<<endl;
            exit(-1);
        }
    } else if (TYPE==Explicit) {
        if (weights.size()!=((N+1)*(N+1)-(N+1))/2) {
            cerr<<"data inconsistency: weights wrong size (Explicit)"<<endl;
            exit(-1);
        }
    } else if (TYPE==FullMatrix) {
        if (weights.size()!=(N+1)*(N+1)) {
            cerr<<"data inconsistency: weights wrong size (FullMatrix)"<<endl;
            exit(-1);
        }
    } else {
        cerr<<"data inconsistency: instance type unspecified"<<endl;
        exit(-1);
    }
    if (probDist!=Louveaux3 && dparams.size()!=N+1) {
        cerr<<"data inconsistency: dparams wrong size"<<endl;
        exit(-1);
    }
    if (probDist==TwoPointDist && (probZero<=5e-2 || probZero>=1-5e-2)) {
        cerr<<"incorrect (?) probZero value (TwoPoint dist.): "<<probZero<<endl;
        exit(-1);
    }
}

int Instance::expDemand(int i) const {
    if (probDist==PoissonDist || probDist==NegBinomialDist
            || probDist==MixedDists || probDist==UniformDist)
        return dparams[i];
    else if (probDist==Louveaux3 || probDist==Louveaux9) {
        return 5;
    } else {
        cerr<<"expDemand: probDist "<<probDist<<" not implemented"<<endl;
        exit(-1);
    }
    return 0;
}

void Instance::genDemandProbs() {
    maxdemands.push_back(-1);                   // depot has no demand
    demandprobs.push_back(vector<double>());
    for (int i=1; i<=N; ++i) {
        vector<double> probs;
        double cum=0;
        int k=0;
        while (cum<MY_ONE) {
            double p=0;
            if (probDist==PoissonDist)
                p=Poisson::pmf(k, dparams.at(i));
            else if (probDist==UniformDist)
                p=1.0/(1+2*dparams.at(i));
            else if (probDist==TwoPointDist) {
                if (k==0)
                    p=probZero;
                else if (k==dparams.at(i))
                    p=1-probZero;
                else
                    p=0;
            } else if (probDist==Louveaux3) {
                if (k==4)
                    p=0.25;
                else if (k==5)      // TODO: typo in the paper formula?
                    p=0.5;
                else if (k==6)
                    p=0.25;
                else
                    p=0.0;
            } else if (probDist==Louveaux9) {
                if (k==1)
                    p=1.0/25;
                else if (k==2)
                    p=2.0/25;
                else if (k==3)
                    p=3.0/25;
                else if (k==4)
                    p=4.0/25;
                else if (k==5)
                    p=5.0/25;
                else if (k==6)
                    p=4.0/25;
                else if (k==7)
                    p=3.0/25;
                else if (k==8)
                    p=2.0/25;
                else if (k==9)
                    p=1.0/25;
            } else if (probDist==NegBinomialDist) {
                p=NegBinomial::pmf(k, dparams.at(i), 0.5);
            } else if (probDist==MixedDists) {
                if (i%3==0) {               // deterministic
                    if (k==dparams[i])
                        p=1.0;
                } else {                    // high variance
                    p=NegBinomial::pmf(k, dparams.at(i), 0.5);
                }
            } else {
                cerr<<"unspecified probability distribution"<<endl;
                exit(-1);
            }
            probs.push_back(p);
            cum+=p;
            k++;
        }
        maxdemands.push_back(probs.size()-1);
        demandprobs.push_back(probs);
    }
}

void Instance::genDistanceMatrix() {
    if (TYPE==XY) {
        for (int i=0; i<xcoords.size(); ++i) {
            vector<double> v;
            for (int j=0; j<ycoords.size(); ++j) {
                double dist=sqrt((xcoords[i]-xcoords[j])
                        *(xcoords[i]-xcoords[j])+(ycoords[i]-ycoords[j])
                        *(ycoords[i]-ycoords[j]));
                if (roundInteger)
                    v.push_back((int)(dist+0.5));
                else
                    v.push_back(dist);
            }
            dmatrix.push_back(v);
        }
    } else if (TYPE==Explicit) {
        dmatrix=vector<vector<double>>(N+1, vector<double>(N+1, 0.0));
        int k=0;
        for (int i=0; i<=N; ++i)
            for (int j=0; j<i; ++j) {
                dmatrix[i][j]=weights[k];
                dmatrix[j][i]=weights[k];
                k++;
            }
        if (k!=weights.size()) {
            cerr<<"genDistanceMatrix(): inconsistency (Explicit)"<<endl;
            exit(-1);
        }
    } else if (TYPE==FullMatrix) {
        dmatrix=vector<vector<double>>(N+1, vector<double>(N+1, 0.0));
        int k=0;
        for (int i=0; i<=N; ++i)
            for (int j=0; j<=N; ++j) {
                if (i!=j)
                    dmatrix[i][j]=weights[k];
                k++;
            }
        if (k!=weights.size()) {
            cerr<<"genDistanceMatrix(): inconsistency (FullMatrix)"<<endl;
            exit(-1);
        }
    } else {
        cerr<<"genDistanceMatrix(): unspecified instance type"<<endl;
        exit(-1);
    }
}

void Instance::initialize() {
    checkData();
    genDistanceMatrix();
    genDemandProbs();
    if (Config::FIX_VEHICLES)
        computeMaxVehicles();
}

void Instance::loadCVRPOptimal(string filename) {
    filename+=".opt";
    ifstream infile(filename);
    if (!infile.is_open()) {
        cout<<"loadCVRPOptimal(): unable to open file "<<filename<<endl;
        cout<<"proceeding anyway..."<<endl;
        return;
    }
    int R;              // number of routes in the optimal solution
    infile>>R;
    for (int i=0; i<R; ++i) {
        vector<int> nodes;
        while (true) {
            int n;
            infile>>n;
            nodes.push_back(n);
            if (n==0 && nodes.size()>1) {
                cvrpRoutes.push_back(nodes);
                break;   // next route
            }
        }
    }
    // validation
    vector<int> visitsNode(N+1, 0);
    for (const auto& r : cvrpRoutes) {
        if (r.size()<=2) {
            cerr<<"loadCVRPOptimal(): invalid route: r.size()=="
                    <<r.size()<<endl;
            cerr<<"CVRP opt file: "<<filename<<endl;
            exit(-1);
        }
        if (r.at(0)!=0 || r.back()!=0) {
            cerr<<"loadCVRPOptimal(): invalid route"<<endl;
            exit(-1);
        }
        for (auto n : r)
            visitsNode.at(n)++;
    }
    for (int i=0; i<visitsNode.size(); ++i)
        if (i!=0 && visitsNode[i]!=1) {
            cerr<<"loadCVRPOptimal(): inconsistency: node "<<i<<" visited "
                    <<visitsNode[i]<<" times"<<endl;
            exit(-1);
        }
}

void Instance::loadEdgeFullMatrixInstance(string filename) {
    TYPE=FullMatrix;
    string filein=filename;
    filein+=".in";
    ifstream infile(filein);
    if (!infile.is_open()) {
        cerr<<"loadEdgeFullMatrixInstance(): unable to open file "
                <<filein<<endl;
        exit(-1);
    }
    infile>>N;
    N--;                    // N does not include the depot
    infile>>Q;
    for (int i=0; i<=N; ++i)
        for (int j=0; j<=N; ++j) {
            int dist;
            infile>>dist;
            weights.push_back(dist);
        }
    dparams=vector<int>(N+1, 0);
    for (int i=0; i<=N; ++i) {
        int n;
        infile>>n>>dparams.at(i);
        if (n!=i+1) {
            cerr<<"loadEdgeFullMatrixInstance("<<filein<<"): inconsistency"
                    <<endl;
            exit(-1);
        }
    }
    loadCVRPOptimal(filename);
}

void Instance::loadEdgeLowerRowInstance(string filename) {
    TYPE=Explicit;
    string filein=filename;
    filein+=".in";
    ifstream infile(filein);
    if (!infile.is_open()) {
        cerr<<"loadEdgeLowerRowInstance(): unable to open file "<<filein<<endl;
        exit(-1);
    }
    infile>>N;
    N--;                    // N does not include the depot
    infile>>Q;
    for (int i=0; i<=N; ++i)
        for (int j=0; j<i; ++j) {
            int dist;
            infile>>dist;
            weights.push_back(dist);
        }
    dparams=vector<int>(N+1, 0);
    for (int i=0; i<=N; ++i) {
        int n;
        infile>>n>>dparams.at(i);
        if (n!=i+1) {
            cerr<<"loadEdgeLowerRowInstance("<<filein<<"): inconsistency"<<endl;
            exit(-1);
        }
    }
    loadCVRPOptimal(filename);
}

void Instance::loadXYInstance(string filename) {
    string filein=filename;
    filein+=".in";
    ifstream infile(filein);
    if (!infile.is_open()) {
        cerr<<"loadXYInstance(): unable to open file "<<filein<<endl;
        exit(-1);
    }
    infile>>N;
    N--;                    // N does not include the depot
    infile>>Q;
    // read x,y coordinates
    xcoords=vector<int>(N+1, 0);
    ycoords=vector<int>(N+1, 0);
    dparams=vector<int>(N+1, 0);
    int n;
    for (int i=0; i<=N; ++i) {
        infile>>n>>xcoords.at(i)>>ycoords.at(i);
        if (n!=i+1) {
            cerr<<"loadXYInstance("<<filein<<"): inconsistency"<<endl;
            exit(-1);
        }
    }
    for (int i=0; i<=N; ++i) {
        infile>>n>>dparams.at(i);
        if (n!=i+1) {
            cerr<<"loadXYInstance("<<filein<<"): inconsistency"<<endl;
            exit(-1);
        }
    }
    loadCVRPOptimal(filename);
}

void Instance::printDistanceMatrix() const {
    for (auto i : dmatrix) {
        for (auto j : i)
            cout<<j<<" ";
        cout<<endl;
    }
}

void Instance::printInfo() const {
    cout<<"Instance info:"<<endl;
    cout<<"ID: "<<ID<<endl;
    cout<<"N: "<<N<<endl;
    cout<<"Q: "<<Q<<endl;
    cout<<"f: "<<f<<endl;
    cout<<"round to integer: "<<roundInteger<<endl;
    if (Config::FIX_VEHICLES)
        cout<<"max vehicles: "<<maxV<<endl;
    else
        cout<<"number of vehicles unfixed"<<endl;
    if (probDist==PoissonDist) {
        cout<<"demand distribution: Poisson"<<endl;
        double avg=0;
        for (auto& d : dparams)
            avg+=d;
        cout<<"expected demand: "<<avg<<endl;
        cout<<"load coeff.: "<<avg/Q<<endl;
    } else if (probDist==UniformDist) {
        cout<<"demand distribution: Discrete Uniform"<<endl;
        int avg=0;
        for (auto& d : dparams)
            avg+=d;
        cout<<"expected demand: "<<avg<<endl;
        cout<<"load coeff.: "<<avg/Q<<endl;
    } else if (probDist==TwoPointDist) {
        cout<<"demand distribution: Two Point"<<endl;
        cout<<"Prob[D=0] = "<<probZero<<endl;
        int max=0;
        for (auto& d : dparams)
            max+=d;
        double avg=(1-probZero)*max;
        cout<<"expected demand: "<<avg<<endl;
        cout<<"load coeff.: "<<avg/Q<<endl;
        cout<<"node data:"<<endl;
        cout<<"<node number> <x coord> <y coord> <dparam>:"<<endl;
        for (int i=0; i<=N; ++i)
            cout<<i<<"\t"<<xcoords[i]<<"\t"<<ycoords[i]<<"\t"<<dparams[i]<<endl;
    } else if (probDist==Louveaux3) {
        cout<<"demand distribution: Louveaux3: [4,5,6] [0.25,0.5,0.25]"<<endl;
        cout<<"expected demand: "<<N*5<<endl;
    } else if (probDist==Louveaux9) {
        cout<<"demand distribution: Louveaux9: [1,2,3,4,5,6,7,8,9] "
                <<"1/25*[1,2,3,4,5,4,3,2,1]"<<endl;
        cout<<"expected demand: "<<N*5<<endl;
    } else if (probDist==NegBinomialDist) {
        cout<<"demand distribution: Negative Binomial"<<endl;
    } else if (probDist==MixedDists) {
        cout<<"demand distribution: mixed"<<endl;
    } else {
        cerr<<"probability distribution unknown"<<endl;
        exit(-1);
    }
}

double Instance::prob(int i, int k) const {
    if (k>maxDemand(i))
        return 0.0;
    return demandprobs.at(i).at(k);
}

