#Compiler
FC = gfortran

#Base option
BFFLAG = -fimplicit-none -std=gnu
BLDFLAG =
BLLFLAG = -llapack -lblas -lpthread
#Run compilation option
RFFLAG = -O2 -march=native -fopenmp
RLDFLAG = -O2 -march=native -fopenmp
#Debug compilation option
DFFLAG = -g -p -pg -W -Wall -pedantic -O0 -g3 -O0 -fno-second-underscore -Wall -Wextra -Wno-compare-reals -fcheck=all -fimplicit-none -std=f2008 -pedantic -ffpe-trap=invalid,zero,overflow -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function -finit-real=snan -finit-logical=false -finit-integer=999999
DLDFLAG = -g -p -pg -O0 -g3 -O0 -fno-second-underscore -Wall -Wextra -Wno-compare-reals -fcheck=all -fimplicit-none -std=f2008 -pedantic -ffpe-trap=zero,overflow -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function

#Default compilation option
FFLAG = $(BFFLAG) $(RFFLAG)
LDFLAG = $(BLDFLAG) $(RLDFLAG)
LLFLAG = $(BLLFLAG)


#Path
SRCDIR = src
OBJDIR = build
BINDIR = bin

#File
#SOURCES := $(wildcard $(SRCDIR)/*.f90)
MOD_SOURCES = constantProgram.f90 constantPhysic.f90 constantSimulation.f90 io.f90 algebricEquations.f90 sCurve.f90 diffSolver.f90 simulationManager.f90
PROG_SOURCES = trouNoir.f90 test/testTemperatureScheme.f90 test/testDensityExplicit.f90 test/testSimulationExplicitScheme.f90 test/testSCurveMesh.f90 test/testingDicho.f90 test/testTransformation.f90 test/testInitialCondition.f90 test/testSCurveProfileTemperature.f90
MOD_SOURCES := $(MOD_SOURCES:%.f90=$(SRCDIR)/%.f90)
PROG_SOURCES := $(PROG_SOURCES:%.f90=$(SRCDIR)/%.f90)
SOURCES := $(MOD_SOURCES) $(PROG_SOURCES)
MOD_OBJECTS := $(MOD_SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
PROG_OBJECTS := $(PROG_SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJECTS := $(SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
MOD := $(SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.mod)



all : $(OBJDIR) $(BINDIR) trouNoir

test : $(OBJDIR) $(BINDIR) testSimulationExplicitScheme testTemperatureScheme testDensityExplicit testingDicho testSCurveMesh testTransformation testInitialCondition testSCurveProfileTemperature

#Construct repertory if doesnt exist
$(OBJDIR):
	mkdir $(OBJDIR)
	mkdir $(OBJDIR)/test

$(BINDIR):
	mkdir $(BINDIR)

#In order to compile with debug option
debug : mrproper
debug : FFLAG = $(BFFLAG) $(DFFLAG)
debug : LDFLAG = $(BSLDFLAG) $(DLDFLAG)
debug : $(OBJDIR) $(BINDIR)
debug : trouNoir test

trouNoir : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/trouNoir.o -o $(BINDIR)/$@ $(LLFLAG)
	cp -f $(BINDIR)/$@ $@

testTemperatureScheme : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testTemperatureScheme.o -o $(BINDIR)/$@ $(LLFLAG)

testDensityExplicit : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testDensityExplicit.o -o $(BINDIR)/$@ $(LLFLAG)

testSimulationExplicitScheme : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testSimulationExplicitScheme.o -o $(BINDIR)/$@ $(LLFLAG)

testSCurveMesh : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testSCurveMesh.o -o $(BINDIR)/$@ $(LLFLAG)

testingDicho : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testingDicho.o -o $(BINDIR)/$@ $(LLFLAG)

testTransformation : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testTransformation.o -o $(BINDIR)/$@ $(LLFLAG)

testInitialCondition : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testInitialCondition.o -o $(BINDIR)/$@ $(LLFLAG)

testSCurveProfileTemperature : $(OBJDIR) $(BINDIR) $(OBJECTS)
	$(FC) $(LDFLAG) $(MOD_OBJECTS) $(OBJDIR)/test/testSCurveProfileTemperature.o -o $(BINDIR)/$@ $(LLFLAG)

#Generic compile command
$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAG) -c $< -o $@ -J$(OBJDIR)

clean :
	rm -rf build/

mrproper : clean
	rm -f trouNoir
	rm -rf bin/
