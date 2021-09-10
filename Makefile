# --- Variable definitions ---
CC = icc
CXX = icc
#CPPFLAGS = -fopenmp -std=c++17 -Wall
#LDFLAGS = -lm -ldl -lgsl -lfftw3 -lgslcblas

CPPFLAGS = -qopenmp -Wall -std=c++17 -fPIC -DMKL_ILP64 -m64 -I"${MKLROOT}/include" 
CPPFLAGS += -I"/home/jalen/.local/include" # If on Boyle
#CPPFLAGS += -g  # activate debugging
#CPPFLAGS += -no-prec-div # approximate division

#CPPFLAGS = -fopenmp -Wall -DMKL_ILP64 -mkl=parallel -I"${MKLROOT}/include" #-qopt-report -qopt-report-phase=openmp
#LDFLAGS = -lm -lgsl -lstdc++ -lfftw3 -Wl,--no-as-needed -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${MKLROOT}/lib/intel64 -ldl
#LDFLAGS = -qopenmp -lfftw3 -liomp5 -L${MKLROOT}/lib/intel64 -lpthread -lm -ldl -lgsl
LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lgsl -lm -ldl
LDFLAGS += -L"/home/jalen/.local/lib" 
#LDFLAGS += -liomp5 -lpthread

OBJFILES = GBPPE.o Materials.o Structure.o createLayers.o Utilities.o
TARGET = test.out
HEADERS = BPPE.h Materials.h physicalConstants.h Structure.h Utilities.h createLayers.h

# --- Rule 1 ---
test: $(TARGET)

$(TARGET): $(OBJFILES) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

rootfind: rootfindtest.cpp
	$(CC) $(CPPFLAGS) -o testrootfind.out rootfindtest.cpp $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
