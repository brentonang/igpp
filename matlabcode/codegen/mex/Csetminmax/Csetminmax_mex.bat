@echo off
set MATLAB=C:\PROGRA~1\MATLAB\MATLAB~1\R2015a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=Csetminmax_mex
set MEX_NAME=Csetminmax_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\MATLAB~1\R2015a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for Csetminmax > Csetminmax_mex.mki
echo COMPILER=%COMPILER%>> Csetminmax_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> Csetminmax_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> Csetminmax_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> Csetminmax_mex.mki
echo LINKER=%LINKER%>> Csetminmax_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> Csetminmax_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> Csetminmax_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> Csetminmax_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> Csetminmax_mex.mki
echo BORLAND=%BORLAND%>> Csetminmax_mex.mki
echo OMPFLAGS= >> Csetminmax_mex.mki
echo OMPLINKFLAGS= >> Csetminmax_mex.mki
echo EMC_COMPILER=lcc64>> Csetminmax_mex.mki
echo EMC_CONFIG=optim>> Csetminmax_mex.mki
"C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\win64\gmake" -B -f Csetminmax_mex.mk
