# --- Variable definitions ---
CC = g++
CXX = g++
CPPFLAGS = -fopenmp -std=c++17 -Wall
LDFLAGS = -lm -ldl -lgsl -lfftw3 -lgslcblas
#CPPFLAGS = -qopenmp -Wall -std=c++17 -fPIC -DMKL_ILP64 -m64 -I"${MKLROOT}/include" #-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3
#CPPFLAGS = -fopenmp -Wall -DMKL_ILP64 -mkl=parallel -I"${MKLROOT}/include" #-qopt-report -qopt-report-phase=openmp
#LDFLAGS = -lm -lgsl -lstdc++ -lfftw3 -Wl,--no-as-needed -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${MKLROOT}/lib/intel64 -ldl
#LDFLAGS = -qopenmp -lfftw3 -liomp5 -L${MKLROOT}/lib/intel64 -lpthread -lm -ldl -lgsl
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lm -ldl
OBJFILES = BPPE.o Materials.o Structure.o
TARGET = test.out
HEADERS = BPPE.h Materials.h physicalConstants.h Structure.h Utilities.h

# --- Rule 1 ---
all: $(TARGET)

$(TARGET): $(OBJFILES) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
