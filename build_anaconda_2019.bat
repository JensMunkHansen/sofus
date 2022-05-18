REM This can be made a lot more elegant - focus is right now to make it work

@echo off

pushd %~dp0

set SWIG_DIR=C:\Program Files\SWIG

REM Only required if documentation is build
set SEDPATH=c:/cygwin64/bin
set GRAPHVIZ_PATH=C:/PFx64/Graphviz2.38/bin
set DOXYGEN_DIR=C:/PFx64/doxygen

REM Only required if tests are build (-DBUILD_FNM_TEST=ON)
set GTEST_ROOT="C:/Program Files/googletest-distribution"

set FFTW_DIR="c:/Program Files/fftw-3.3.5"

REM Visual Studio environment
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64

REM Python virtual or base environment
call "c:\ProgramData\Anaconda3\Scripts\activate.bat"

mkdir build

set PATH="c:/PFx64/CMake/bin";%SWIG_DIR%;%SEDPATH%;%GRAPHVIZ_PATH%;%PATH%

REM Create solution
cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 16 2019" -A "x64" -DBUILD_SPS_TEST=OFF -DBUILD_DOCUMENTATION=OFF -DBUILD_SWIG_DOCUMENTATION=OFF -DBUILD_FNM_TEST=OFF -DBUILD_GL_TEST=OFF -DDOXYGEN_EXECUTABLE=%DOXYGEN_DIR%\bin\doxygen.exe -DSOFUS_DEBUG_USING_PYTHON_RELEASE_RUNTIME=ON

cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

cd %~dp0
