SYSTEM = x86-64_osx
LIBFORMAT = static_pic

CPLEXDIR = /Applications/CPLEX_Studio1210/cplex
CONCERTDIR = /Applications/CPLEX_Studio1210/concert

CCC = clang++ -O2 -Wall -std=c++11
#CCOPT = -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD -stdlib=libc++ -c
CCOPT = -m64 -fPIC -fexceptions -DIL_STD -stdlib=libc++ -c

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -L./cvrpsep
CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit -stdlib=libc++ -lcvrpsep

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

SOURCES=Arc.cpp Binomial.cpp Config.cpp Data.cpp Instance.cpp IntegerSolver.cpp Knapsack.cpp Label.cpp LinearSolver.cpp NegBinomial.cpp Poisson.cpp PolicySolverI.cpp Pricing.cpp RoundedCapacityCut.cpp Route.cpp Separation3SRC.cpp SeparationRCC.cpp Solution.cpp Solver.cpp Stats.cpp SubsetRowCut.cpp SwitchPolicy.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CCC) $(CCLNDIRS) $(CCLNFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CCC) $(CCFLAGS) $< -o $@

clean:
	/bin/rm -rf *.o

