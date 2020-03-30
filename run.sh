export MATLAB="/Applications/MATLAB_R2019b.app/"
export DYLD_LIBRARY_PATH=$MATLAB/bin/maci64/:$MATLAB/sys/os/maci64/:$DYLD_LIBRARY_PATH

# add matlab engine
export PATH=/Applications/MATLAB_R2019b.app/bin:$PATH

make -f Makefile

./run	
