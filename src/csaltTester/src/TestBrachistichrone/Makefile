# Project:        Collocation
# Created by:     Jeremy Knittel
# Date:           2016.11.03


# Define macros for locations
BASE_ROOT = /Users/granthec/Source/CSALT/src/csalt
TEST_ROOT = /Users/granthec/Source/CSALT/src/csaltTester
GMAT_Utilities = /Users/granthec/Source/CSALT/src/gmatutil
BOOST_ROOT = /Users/granthec/Source/CSALT/depends/Boost_1_71_0
SNOPT_ROOT = /Users/granthec/Source/CSALT/depends/snopt7
BINDIR = $(TEST_ROOT)/bin
EXECUTABLE = $(BINDIR)/TestBrachistichrone

# TEST_FLAGS = -DUSE_64_BIT_LONGS -D__CONSOLE_APP__
TEST_FLAGS = -DUSE_64_BIT_LONGS -D__DEBUG__ -D__CONSOLE_APP__

#OPTIMIZATIONS = -O3 -funroll-loops
OPTIMIZATIONS = 

# Define macros for the needed includes
HEADERS =   -I. \
            -I$(TEST_ROOT)/src/TestBrachistichrone \
            -I$(TEST_ROOT)/src/HelperClasses \
            -I$(BASE_ROOT)/src/userfunutils \
            -I$(BASE_ROOT)/src/collutils \
            -I$(BASE_ROOT)/src/executive \
            -I$(BASE_ROOT)/src/util \
            -I$(BASE_ROOT)/src/include \
            -I$(GMAT_Utilities)/include \
            -I$(GMAT_Utilities)/util \
            -I$(GMAT_Utilities)/util/datawriter \
            -I$(GMAT_Utilities)/util/interpolator \
            -I$(GMAT_Utilities)/util/matrixoperations \
            -I$(BOOST_ROOT) \
            -I$(SNOPT_ROOT)/include \
            -I$(SNOPT_ROOT)/src \
            
DRIVER = TestBrachistichrone.o
            
OBJECTS = \
    $(TEST_ROOT)/src/HelperClasses/BrachistichronePathObject.o \
    $(TEST_ROOT)/src/HelperClasses/BrachistichronePointObject.o \
    $(BASE_ROOT)/src/util/LowThrustException.o \
    $(BASE_ROOT)/src/util/RadauMathUtil.o \
    $(BASE_ROOT)/src/util/BaryLagrangeInterpolator.o \
    $(BASE_ROOT)/src/util/SparseMatrixUtil.o \
	$(BASE_ROOT)/src/util/SNOPTFunctionWrapper.o \
    $(BASE_ROOT)/src/executive/Phase.o \
    $(BASE_ROOT)/src/executive/RadauPhase.o \
    $(BASE_ROOT)/src/executive/OptimizerBase.o \
    $(BASE_ROOT)/src/executive/Trajectory.o \
    $(BASE_ROOT)/src/collutils/DecisionVector.o \
    $(BASE_ROOT)/src/collutils/DecVecTypeBetts.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtilRadau.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtil.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtil_AlgPath.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtil_Coll.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtil_MultiPoint.o \
    $(BASE_ROOT)/src/collutils/NLPFuncUtil_Path.o \
    $(BASE_ROOT)/src/collutils/NLPFunctionData.o \
    $(BASE_ROOT)/src/collutils/ProblemCharacteristics.o \
    $(BASE_ROOT)/src/userfunutils/BoundData.o \
    $(BASE_ROOT)/src/userfunutils/JacobianData.o \
    $(BASE_ROOT)/src/userfunutils/UserFunction.o \
    $(BASE_ROOT)/src/userfunutils/UserFunctionManager.o \
    $(BASE_ROOT)/src/userfunutils/UserPathFunctionManager.o \
    $(BASE_ROOT)/src/userfunutils/UserPointFunctionManager.o \
    $(BASE_ROOT)/src/userfunutils/UserFunctionProperties.o \
    $(BASE_ROOT)/src/userfunutils/UserPathFunction.o \
    $(BASE_ROOT)/src/userfunutils/UserPointFunction.o \
    $(BASE_ROOT)/src/userfunutils/FunctionOutputData.o \
    $(BASE_ROOT)/src/userfunutils/FunctionInputData.o \
    $(BASE_ROOT)/src/userfunutils/FunctionContainer.o \
    $(BASE_ROOT)/src/userfunutils/PathFunctionContainer.o \
    $(BASE_ROOT)/src/userfunutils/PointFunctionContainer.o \
    $(GMAT_Utilities)/util/BaseException.o \
    $(GMAT_Utilities)/util/FileManager.o \
    $(GMAT_Utilities)/util/FileUtil.o \
    $(GMAT_Utilities)/util/GmatGlobal.o \
    $(GMAT_Utilities)/util/Linear.o \
    $(GMAT_Utilities)/util/MessageInterface.o \
    $(GMAT_Utilities)/util/MessageReceiver.o \
    $(GMAT_Utilities)/util/RandomNumber.o \
    $(GMAT_Utilities)/util/RealUtilities.o \
    $(GMAT_Utilities)/util/Rmatrix.o \
    $(GMAT_Utilities)/util/Rmatrix33.o \
    $(GMAT_Utilities)/util/Rvector.o \
    $(GMAT_Utilities)/util/Rvector3.o \
    $(GMAT_Utilities)/util/StringTokenizer.o \
    $(GMAT_Utilities)/util/StringUtil.o \
    $(GMAT_Utilities)/util/TimeTypes.o

###    $(BASE_ROOT)/src/userfunutils/DummyPathFunction.o \

CPPFLAGS = $(HEADERS) $(OPTIMIZATIONS) $(TESTFLAGS) -Wall -std=c++17
CPP = g++

# Link specific flags
LINK_FLAGS =  -headerpad_max_install_names \
               -lm -v \
			   -L$(SNOPT_ROOT)/lib/.libs/ \
			   -lsnopt7_cpp \
			   -lsnopt7

# Object modules for library
$(DRIVER): %.o: %.cpp
	$(CPP) $(CPPFLAGS) $(HEADERS) -c -o $@ $<

$(OBJECTS): %.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(HEADERS) -c -o $@ $<

#
# Default suffix rule
.SUFFIXES: .o .cpp

# Target definitions
#
all: $(EXECUTABLE)

$(EXECUTABLE): $(DRIVER) $(OBJECTS)
	@if [ ! -x $(BINDIR) ]; then mkdir $(BINDIR); fi 
	$(CPP) $(CPPFLAGS) $(LINK_FLAGS) $(DRIVER) $(OBJECTS) -o $@
	@chmod -R 755 $(BINDIR)

clean:
	rm -f $(DRIVER) $(OBJECTS) core

bare: clean
	rm -f $(EXECUTABLE) 
