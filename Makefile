CC = g++
CFLAGS = -std=c++11
INCLUDE = -I/Applications/MATLAB_R2019b.app/extern/include/ -I/Applications/MATLAB_R2019b.app/simulink/include
LIB1 = /Applications/MATLAB_R2019b.app/bin/maci64
LIB2 = /Applications/MATLAB_R2019b.app/extern/bin/maci64

object_file = run.o

run	:	$(object_file)
	$(CC) $(CFLAGS) $(object_file) -o run -I$(INCLUDE) -L$(LIB1) -L$(LIB2) -lmx -lmex -lmat -leng -lMatlabDataArray -lMatlabEngine

run.o	:	run.cpp MatLib.h
	$(CC) -c run.cpp -o run.o $(INCLUDE)

clean	:
	rm $(object_file) run
