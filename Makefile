INCPATH = ./include
SRC = lib/blas.cpp lib/dirac_op.cpp lib/hmc.cpp lib/inverters.cpp \
	lib/io.cpp lib/measurements.cpp lib/utils.cpp

CXX = g++
CXXFLAGS = -O3 -Xpreprocessor -fopenmp -std=c++0x
#CXXFLAGS = -static-libstdc++ -O3 -fopenmp -std=c++0x

all: schwinger2peD testDirac testHMC

schwinger2peD: schwinger2peD.cpp $(SRC) $(INCPATH)/* Makefile
	@echo Compiling schwinger2peD...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(SRC) schwinger2peD.cpp -o $@ -lomp

testDirac: testDirac.cpp $(SRC) $(INCPATH)/* Makefile
	@echo Compiling testDirac...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(SRC) testDirac.cpp -o $@ -lomp

testHMC: testHMC.cpp $(SRC) $(INCPATH)/* Makefile
	@echo Compiling testHMC...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(SRC) testHMC.cpp -o $@ -lomp
