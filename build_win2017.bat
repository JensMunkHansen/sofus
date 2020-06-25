@echo off

pushd %~dp0

call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

set FFTWDIR="C:\Program Files\fftw-3.3.5"

mkdir build

set SWIG_DIR=C:\Program Files\swig

REM Ensure the right version of CMake is used
set PATH="C:\Program Files\CMake\bin";%SWIG_DIR%;%PATH%

REM cmake -H%~dp0\deps -B%~dp0\build\deps-prefix -G "Visual Studio 15 2017 Win64"

REM cd build\deps-prefix
REM msbuild fnm_deps.sln /p:Configuration=Release

cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 15 2017 Win64" -DBUILD_SPS_TEST=ON -DBUILD_FNM_TEST=ON -DBUILD_GL_TEST=ON

REM Create solution
cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

REM  /verbosity:diagnostic
cd ..

cd %~dp0
