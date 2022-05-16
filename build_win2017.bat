@echo off

pushd %~dp0

set SWIG_DIR=C:\Program Files\SWIG

REM Only required if documentation is build
set SEDPATH=c:/cygwin64/bin
set GRAPHVIZ_PATH=C:\Program Files (x86)\Graphviz2.38\bin
set DOXYGEN_DIR=C:/PFx64/doxygen
set CONDAPATH=C:\ProgramData\Anaconda3

REM Only required if tests are build (-DBUILD_FNM_TEST=ON)
set GTEST_ROOT="C:/Program Files/googletest-distribution"

call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

REM Python virtual or base environment
call "c:\ProgramData\Anaconda3\Scripts\activate.bat"

mkdir build

REM Ensure the right version of CMake is used
set PATH="C:\Program Files\CMake\bin";%SWIG_DIR%;%SEDPATH%;%GRAPHVIZ_PATH%;%PATH%

REM cmake -H%~dp0\deps -B%~dp0\build\deps-prefix -G "Visual Studio 15 2017 Win64"

REM cd build\deps-prefix
REM msbuild fnm_deps.sln /p:Configuration=Release

cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 15 2017 Win64" -DBUILD_SPS_TEST=OFF -DBUILD_DOCUMENTATION=OFF -DBUILD_FNM_TEST=OFF -DBUILD_GL_TEST=OFF -DPYTHON_INCLUDE_DIR=%CONDAPATH%\include -DPYTHON_LIBRARY=%CONDAPATH%\libs\python37.lib -DPYTHON_EXECUTABLE=%CONDAPATH%\python.exe

REM Create solution
cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

REM  /verbosity:diagnostic
cd ..

cd %~dp0
