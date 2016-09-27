# sofus
SIMD optimized fast ultrasound simulator

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


