MATLAB="/afs/ece.cmu.edu/support/matlab/2013a"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/afs/ece/usr/gustavor/.matlab/R2013a"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for GPReduce" > GPReduce_mex.mki
echo "CC=$CC" >> GPReduce_mex.mki
echo "CFLAGS=$CFLAGS" >> GPReduce_mex.mki
echo "CLIBS=$CLIBS" >> GPReduce_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> GPReduce_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> GPReduce_mex.mki
echo "CXX=$CXX" >> GPReduce_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> GPReduce_mex.mki
echo "CXXLIBS=$CXXLIBS" >> GPReduce_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> GPReduce_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> GPReduce_mex.mki
echo "LD=$LD" >> GPReduce_mex.mki
echo "LDFLAGS=$LDFLAGS" >> GPReduce_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> GPReduce_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> GPReduce_mex.mki
echo "Arch=$Arch" >> GPReduce_mex.mki
echo OMPFLAGS= >> GPReduce_mex.mki
echo OMPLINKFLAGS= >> GPReduce_mex.mki
echo "EMC_COMPILER=" >> GPReduce_mex.mki
echo "EMC_CONFIG=optim" >> GPReduce_mex.mki
"/afs/ece.cmu.edu/support/matlab/2013a/bin/glnxa64/gmake" -B -f GPReduce_mex.mk
