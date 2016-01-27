CC=g++
CFLAGS=-g 

all: fdtdWave1DModel  

fdtdWave1DModel: fdtdWave1DModel.cpp fdtdAcoustic_aux.cpp fdtdViscosity_aux.cpp  fdtdBoundaries_aux.cpp fdtdWave_aux.cpp Array_aux.cpp IO_utils.cpp
	${CC} -o fdtdWave1DModel ${CFLAGS} fdtdWave1DModel.cpp fdtdAcoustic_aux.cpp fdtdViscosity_aux.cpp  fdtdBoundaries_aux.cpp fdtdWave_aux.cpp Array_aux.cpp IO_utils.cpp

	rm -rf *dSYM
