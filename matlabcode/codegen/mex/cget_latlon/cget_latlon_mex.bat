@echo off
set MATLAB=C:\PROGRA~1\MATLAB\MATLAB~1\R2015a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=cget_latlon_mex
set MEX_NAME=cget_latlon_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\MATLAB~1\R2015a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for cget_latlon > cget_latlon_mex.mki
echo COMPILER=%COMPILER%>> cget_latlon_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> cget_latlon_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> cget_latlon_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> cget_latlon_mex.mki
echo LINKER=%LINKER%>> cget_latlon_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> cget_latlon_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> cget_latlon_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> cget_latlon_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> cget_latlon_mex.mki
echo BORLAND=%BORLAND%>> cget_latlon_mex.mki
echo OMPFLAGS= >> cget_latlon_mex.mki
echo OMPLINKFLAGS= >> cget_latlon_mex.mki
echo EMC_COMPILER=lcc64>> cget_latlon_mex.mki
echo EMC_CONFIG=optim>> cget_latlon_mex.mki
"C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\win64\gmake" -B -f cget_latlon_mex.mk
