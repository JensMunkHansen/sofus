sofus
=====

SIMD optimized fast ultrasound simulator. Only continuous-wave pressure
is optimized using SIMD in this repository

## Installation
At the moment, a reference implementation of the fast nearfield method
(FNM) is available for download. The code can be built by issuing

    cd sofus
    mkdir build
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=Release
	make

An example of the continuous-wave (CW) pressure in front of an
128-element array is given in ```test_cw_pressure.py```

    cd ../sofus
	ipython
	run test_cw_pressure.py

Further examples with transients and pulsed-wave fields can be found
in the examples folder.

### FFTW

An implementation of the angular spectrum approach (ASA) for field propagation
is included and this requires FFTW.

#### Windows

For windows, binaries can be fetched from
http://www.fftw.org/install/windows.html . You need to use the lib.exe
tool for generating lib files as described on the webpage. I have put the binaries in C:\Program Files\fftw-3.3.5.


## Documentation
The documentation can be found as the io page. [Doxygen documentation](http://jensmunkhansen.github.io/sofus "SOFUS")
 
