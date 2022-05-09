REM This can be made a lot more elegant - focus is right to make it work

@echo off

pushd %~dp0


set SEDPATH=c:/cygwin64/bin
set GRAPHVIZ_PATH=C:\Program Files (x86)\Graphviz2.38\bin
set SWIG_DIR=C:\Program Files\SWIG
REM REM set SWIG_DIR=C:\Program Files\swig
set DOXYGEN_DIR=C:/PFx64/doxygen

set(GTEST_ROOT="C:/Program Files/googletest-distribution")

set FFTWDIR="C:\Program Files\fftw-3.3.5"

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64

call "c:\cygwin64\home\jem\Environments\PyPI38\Scripts\activate.bat"

mkdir build

set PATH="C:\Program Files\CMake\bin";%SWIG_DIR%;%SEDPATH%;%GRAPHVIZ_PATH%;%PATH%

REM ;%DOXYGEN_DIR%;%GRAPHVIZ_PATH%;%SWIG_DIR%;%PATH%;%SEDPATH%

cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 16 2019" -A "x64" -DBUILD_SPS_TEST=OFF -DBUILD_DOCUMENTATION=ON -DBUILD_SWIG_DOCUMENTATION=ON -DBUILD_FNM_TEST=ON -DBUILD_GL_TEST=ON -DDOXYGEN_EXECUTABLE=%DOXYGEN_DIR%/bin/doxygen.exe -DPYTHON_EXECUTABLE="C:\cygwin64\home\jem\Environments\PyPI38\Scripts\python.exe"

REM Create solution
cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

REM  /verbosity:diagnostic
cd ..

cd %~dp0
