REM This can be made a lot more elegant - focus is right to make it work

@echo off

pushd %~dp0

set CONDAPATH=C:\ProgramData\Anaconda2\envs\Python37

REM set QTDIR=C:\Qt\Qt5.12.9\5.12.9\msvc2017_64
REM set PYPIPATH=c:/Users/jem/Environments/fis
REM set PYPIPATH=c:/ProgramData/Anaconda3
REM REM set PYPIPATH=c:/Users/jem/AppData/Local/Programs/Python/Python36
REM set PYVER=37
REM set SEDPATH=c:/cygwin64/bin
REM set GRAPHVIZ_PATH=C:\Program Files (x86)\Graphviz\bin
REM set SWIG_DIR=C:\Program Files\SWIG4.0
REM REM set SWIG_DIR=C:\Program Files\swig
REM set DOXYGEN_DIR=C:\Program Files\doxygen\bin

set(GTEST_ROOT="C:/Program Files/googletest-distribution")

set FFTWDIR="C:\Program Files\fftw-3.3.5"

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64


mkdir build

set PATH="C:\Program Files\CMake\bin";%SWIG_DIR%;%CONDAPATH%;%PATH%

REM ;%DOXYGEN_DIR%;%GRAPHVIZ_PATH%;%SWIG_DIR%;%PATH%;%SEDPATH%

cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 16 2019" -A "x64" -DBUILD_SPS_TEST=ON -DBUILD_DOCUMENTATION=OFF -DBUILD_FNM_TEST=ON -DBUILD_GL_TEST=ON -DPYTHON_INCLUDE_DIR=%CONDAPATH%\include -DPYTHON_LIBRARY=%CONDAPATH%\libs\python37.lib -DPYTHON_EXECUTABLE=%CONDAPATH%\python.exe

REM Create solution
cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

REM  /verbosity:diagnostic
cd ..

cd %~dp0
