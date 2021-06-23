# --- Variable definitions ---
CC = icpc
CPPFLAGS = -Wall -std=c++17 -fPIC -fopenmp -DMKL_ILP64  -m64  -I"${MKLROOT}/include" #-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3
LDFLAGS = -lm -lgsl -lstdc++ -lfftw3 -Wl,--no-as-needed -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${MKLROOT}/lib/intel64 -ldl
OBJFILES = BPPE.o Materials.o Structure.o
TARGET = test.out
HEADERS = BPPE.h Materials.h physicalConstants.h Structure.h Utilities.h

# --- Rule 1 ---
all: $(TARGET)

$(TARGET): $(OBJFILES) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
