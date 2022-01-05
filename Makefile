# --- Variable definitions ---
CC = icc
CXX = icc

# --- General flags for compiler and linker ---
CPPFLAGS = -qopenmp -Wall -I"${MKLROOT}/include" -fPIC -DMKL_ILP64 -m64 
CXXFLAGS = -std=c++17 
#CPPFLAGS += -I"/home/jalen/.local/include" # If on Boyle
#CPPFLAGS += -g  # activate debugging
#CPPFLAGS += -no-prec-div # approximate division

LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lgsl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lm -ldl
#LDFLAGS += -L"/home/jalen/.local/lib"  # If on Boyle

# --- Custom variables ---
OBJFILES = GBPPE.o Materials.o Structure.o createLayers.o Utilities.o func.o output.o gslStructs.o helperFuncs.o
TARGET = test.out
HEADERS = BPPE.h Materials.h physicalConstants.h Structure.h Utilities.h createLayers.h gslStructs.h helperFuncs.h

# --- Variables for compiling custom GSL library ---
GSL_MULTIROOTOBJS = gsl/multiroot/convergence.o gsl/multiroot/dogleg.o gsl/multiroot/enorm.o gsl/multiroot/fsolver.o gsl/multiroot/hybrid.o gsl/multiroot/broyden.o gsl/multiroot/dnewton.o gsl/multiroot/fdjac.o
GSLHEADERS = gsl/config.h gsl/gsl_math.h gsl/gsl_types.h gsl/multiroot/gsl_multiroots.h

GSLOBJS = $(GSL_MULTIROOTOBJS)

MYGSLLIBS = libmultroot.a
MYLIBFLAGS = -L. -lmultroot


# --- Rule 1 ---
$(TARGET): $(OBJFILES) $(MYGSLLIBS) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(GSLOBJS) $(MYGSLLIBS) $(LDFLAGS)

testWOmygsl: $(OBJFILES) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

mygsl: $(GSLOBJS) $(GSLHEADERS)
	$(CC) -o mygslobj.o $(GSLOBJS) $(LDFLAGS)

libmultroot.a: $(GSL_MULTIROOTOBJS) $(GSLHEADERS)
	ar $(ARFLAGS) $@ $^

rootfind: rootfindtest.cpp
	$(CC) $(CPPFLAGS) -o testrootfind.out rootfindtest.cpp $(LDFLAGS)


# --- Rules for cleaning ---
cleanall:
	rm -f $(OBJFILES) $(TARGET) $(GSLOBJS) $(MYGSLLIBS) *~

clean:
	rm -f $(OBJFILES) $(TARGET) *~

cleangsl:
	rm -f $(GSLOBJS) $(MYGSLLIBS) *~