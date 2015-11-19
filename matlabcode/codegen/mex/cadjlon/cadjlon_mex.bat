@echo off
set MATLAB=C:\PROGRA~1\MATLAB\MATLAB~1\R2015a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=cadjlon_mex
set MEX_NAME=cadjlon_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\MATLAB~1\R2015a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for cadjlon > cadjlon_mex.mki
echo COMPILER=%COMPILER%>> cadjlon_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> cadjlon_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> cadjlon_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> cadjlon_mex.mki
echo LINKER=%LINKER%>> cadjlon_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> cadjlon_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> cadjlon_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> cadjlon_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> cadjlon_mex.mki
echo BORLAND=%BORLAND%>> cadjlon_mex.mki
echo OMPFLAGS= >> cadjlon_mex.mki
echo OMPLINKFLAGS= >> cadjlon_mex.mki
echo EMC_COMPILER=lcc64>> cadjlon_mex.mki
echo EMC_CONFIG=optim>> cadjlon_mex.mki
"C:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\win64\gmake" -B -f cadjlon_mex.mk
