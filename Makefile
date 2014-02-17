# Alignment of Supercooled, Symmetric Top Molecules in FDTD fields
# Josh Szekely, January 2014
# Accepts FDTD file as input
# Loads data, finds rotational PES, trajectory of molecules
# Assumes the FDTD calculation has already been calibrated
# Makefile

CC = icpc

MKLROOT = /opt/intel/composer_xe_2013_sp1/mkl

CFLAGS = -g -debug -O2 -openmp -I$(MKLROOT)/include -mkl

GSL_INC = -I/Users/joshuaszekely/Desktop/NUResearch/Codes/gsl-1.15
BOOST_INC = -I/opt/local/include/

#FFTW_DIR = /Users/joshuaszekely/Desktop/NUResearch/Codes/fftw-3.3.3
#FFTW_LIB = -L$(FFTW_DIR)
#FFTW_INC = -I$(FFTW_DIR)/include

OBJ  =   obj/main.o  obj/array_structs.o  obj/numerics.o obj/quantum_dynamics.o

LIBS = -lm -lgsl
#LIBS =   -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -lgsl

HEADS =  include/array_structs.h include/numerics.h include/quantum_dynamics.h
BIN  =   MolecAlignment

RM = rm  -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(GSL_INC) -o $(BIN) -I./include $(OBJ) $(BOOST_INC) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) $(GSL_INC) -c  src/main.cpp -o obj/main.o $(BOOST_INC) -I./include 

obj/array_structs.o: src/array_structs.cpp
	$(CC) $(CFLAGS)  -c src/array_structs.cpp  -o obj/array_structs.o  -I./include 

obj/numerics.o: src/numerics.cpp
	$(CC) $(CFLAGS) $(GSL_INC) -c src/numerics.cpp  -o obj/numerics.o  -I./include 

obj/quantum_dynamics.o: src/quantum_dynamics.cpp 
	$(CC) $(CFLAGS) -c src/quantum_dynamics.cpp -o obj/quantum_dynamics.o -I./include $(BOOST_INC)
