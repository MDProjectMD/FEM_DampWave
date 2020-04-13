CC = g++
CFLAGS = -std=c++11

INCLUDE = /Applications/MATLAB_R2019b.app/extern/include/ -I/Applications/MATLAB_R2019b.app/simulink/include
LIB1 = /Applications/MATLAB_R2019b.app/bin/maci64
LIB2 = /Applications/MATLAB_R2019b.app/extern/bin/maci64

INCLUDE_EIGEN = /Users/dengxiaohui/Documents/Eigen
INCLUDE_SPECTRA = /Users/dengxiaohui/Documents/Spectra/include/

object_file = run.o Define.o Random.o SystemDefine.o Stochastic.o TimeDependentSolver.o MatEng.o

#create_velocity	:	$(object_file)
#	$(CC) $(CFLAGS) $(object_file) -o create_velocity -I$(INCLUDE) -L$(LIB1) -L$(LIB2) -lmx -lmex -lmat -leng -lMatlabDataArray -lMatlabEngine
run	:	$(object_file)
	$(CC) $(CFLAGS) $(object_file) -o run -I$(INCLUDE) -L$(LIB1) -L$(LIB2) -lmx -lmex -lmat -leng -lMatlabDataArray -lMatlabEngine

run.o	:	run.cpp  Stochastic.h Define.h MatEng.h Matrix.h
	$(CC) -c run.cpp -o run.o -I$(INCLUDE) -I$(INCLUDE_EIGEN)
#create_velocity.o	:	create_velocity.cpp
#	$(CC) -c create_velocity.cpp -o create_velocity.o

Define.o	:	Define.h Define.cpp
	$(CC) -c Define.cpp -o Define.o

Random.o	:	Random.h Random.cpp
	$(CC) -c Random.cpp -o Random.o

SystemDefine.o	:	SystemDefine.h SystemDefine.cpp
	$(CC) -c SystemDefine.cpp -o SystemDefine.o -I$(INCLUDE_EIGEN)

Stochastic.o	:	Stochastic.h Stochastic.cpp
	$(CC) -c Stochastic.cpp -o Stochastic.o

TimeDependentSolver.o	:	TimeDependentSolver.cpp TimeDependentSolver.h
	$(CC) -c TimeDependentSolver.cpp -o TimeDependentSolver.o -I$(INCLUDE_EIGEN) -I$(INCLUDE)

MatEng.o	:	MatEng.h MatEng.cpp
	$(CC) -c MatEng.cpp -o MatEng.o -I$(INCLUDE) -I$(INCLUDE_EIGEN)

clean	:
	rm $(object_file) run
