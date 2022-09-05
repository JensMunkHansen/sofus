REM This can be made a lot more elegant - focus is right now to make it work

@echo off

pushd %~dp0

REM Required for documentation and inline documentation in Python
set SEDPATH=c:/cygwin64/bin
set GRAPHVIZ_PATH=C:\Program Files (x86)\Graphviz2.38\bin
set DOXYGEN_DIR=C:/PFx64/doxygen

set SWIG_DIR=C:\Program Files\swig

set GTEST_ROOT="C:/Program Files/googletest-distribution"

set FFTW_DIR="c:/Program Files/fftw-3.3.5"

REM Visual Studio environment
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" x64

REM Python virtual or base environment
IF EXIST "c:\cygwin64\home\jem\Environments\PyPI38\Scripts\activate.bat" (
  call "c:\cygwin64\home\jem\Environments\PyPI38\Scripts\activate.bat"
) ELSE IF EXIST "C:\Users\Jens Munk Hansen\Environments\Python310\Scripts\activate.bat" (
  call "C:\Users\Jens Munk Hansen\Environments\Python310\Scripts\activate.bat"
) ELSE (
  call "c:\Users\jem\Environments\PyPI38\Scripts\activate.bat"
)   


mkdir build

set PATH="C:\Program Files\CMake\bin";%SWIG_DIR%;%SEDPATH%;%GRAPHVIZ_PATH%;%PATH%;%FFTW_DIR%

REM Create solution
cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 16 2019" -A "x64" -DBUILD_SPS_TEST=OFF -DBUILD_DOCUMENTATION=OFF -DBUILD_SWIG_DOCUMENTATION=OFF -DBUILD_FNM_TEST=OFF -DBUILD_GL_TEST=OFF -DDOXYGEN_EXECUTABLE=%DOXYGEN_DIR%/bin/doxygen.exe

cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

cd %~dp0
