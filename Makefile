INCPATH = ./include
SRC = lib/blas.cpp lib/dirac_op.cpp lib/hmc.cpp lib/inverters.cpp \
	lib/io.cpp lib/measurements.cpp lib/utils.cpp
OBJ = $(SRC:cpp=o)
ALLEXEC = schwinger2peD testDirac testHMC testIO testPlaq testWilsonFlow

CXX = g++
CXXFLAGS = -static-libstdc++ -O3 -std=c++0x

all: $(ALLEXEC)

$(OBJ): $(SRC) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) -c $(@:o=cpp) -o $@

schwinger2peD: schwinger2peD.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

testDirac: testDirac.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

testHMC: testHMC.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

testIO: testIO.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

testPlaq: testPlaq.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

testWilsonFlow: testWilsonFlow.cpp $(OBJ) $(INCPATH)/* Makefile
	@echo Compiling $@...
	$(CXX) $(CXXFLAGS) -I $(INCPATH) $(OBJ) $@.cpp -o $@

clean:
	rm $(OBJ)
	rm $(ALLEXEC)
