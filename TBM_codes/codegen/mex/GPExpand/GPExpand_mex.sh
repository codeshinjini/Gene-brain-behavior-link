MATLAB="/afs/ece.cmu.edu/support/matlab/2013a"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/afs/ece/usr/gustavor/.matlab/R2013a"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for GPExpand" > GPExpand_mex.mki
echo "CC=$CC" >> GPExpand_mex.mki
echo "CFLAGS=$CFLAGS" >> GPExpand_mex.mki
echo "CLIBS=$CLIBS" >> GPExpand_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> GPExpand_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> GPExpand_mex.mki
echo "CXX=$CXX" >> GPExpand_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> GPExpand_mex.mki
echo "CXXLIBS=$CXXLIBS" >> GPExpand_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> GPExpand_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> GPExpand_mex.mki
echo "LD=$LD" >> GPExpand_mex.mki
echo "LDFLAGS=$LDFLAGS" >> GPExpand_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> GPExpand_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> GPExpand_mex.mki
echo "Arch=$Arch" >> GPExpand_mex.mki
echo OMPFLAGS= >> GPExpand_mex.mki
echo OMPLINKFLAGS= >> GPExpand_mex.mki
echo "EMC_COMPILER=" >> GPExpand_mex.mki
echo "EMC_CONFIG=optim" >> GPExpand_mex.mki
"/afs/ece.cmu.edu/support/matlab/2013a/bin/glnxa64/gmake" -B -f GPExpand_mex.mk
