@echo off

pushd %~dp0

call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

set BuildGui=0
if "%~1"=="NOGUI" goto :NoGUI
set BuildGui=0
:NoGUI

set CMAKE_PREFIX_PATH=C:\Qt\Qt5.9.1\5.9.1\msvc2017_64

REM Variable used in this script - name is arbitrary
set QTDIR=C:\Qt\Qt5.9.1\5.9.1\msvc2017_64

mkdir build
mkdir build\deps-prefix

REM Ensure the right version of CMake is used
set PATH="C:\Program Files\CMake\bin";%PATH%

cmake -H%~dp0\deps -B%~dp0\build\deps-prefix -G "Visual Studio 15 2017 Win64"

cd build\deps-prefix
msbuild fnm_deps.sln /p:Configuration=Release

if defined BuildGUI (
    cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 15 2017 Win64" -DBUILD_SPS_TEST=ON -DBUILD_FNM_TEST=ON -DBUILD_GL_TEST=ON -DUSE_ProgressBar=ON -DBUILD_SOFUS_UI=ON
)
if not defined BuildGUI (
    cmake -H%~dp0 -B%~dp0\build -G "Visual Studio 15 2017 Win64" -DBUILD_SPS_TEST=ON -DBUILD_FNM_TEST=ON -DBUILD_GL_TEST=ON -DUSE_ProgressBar=ON -DBUILD_SOFUS_UI=OFF
)

REM Create solution
cd %~dp0\build

REM Build solution
msbuild fnm.sln /p:Configuration=Release

REM  /verbosity:diagnostic
cd ..

REM Go to output directory
cd %~dp0\build\gui\Release

REM Deploy all necessary DLL's to output directory.
if defined BuildGUI (
  set CURRENTDRIVE=%CD:~0,2%
  %QTDIR%\bin\qtenv2.bat
  call %CURRENTDRIVE%
  cd %~dp0\build\gui\Release
  %QTDIR%\bin\windeployqt.exe SofusUI.exe
)

REM Add FFTW
set FFTWDIR="C:\Program Files\fftw-3.3.4"
cp %FFTWDIR%\libfftw3-3.dll %~dp0\build\gui\Release
cp %FFTWDIR%\libfftw3f-3.dll %~dp0\build\gui\Release

cd %~dp0

