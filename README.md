
<h1>FiDiTi - FDTD Library</h1>

Header-only C++ library, implementing the Finite-Difference Time-Domain Method (FDTD)

Requires C++17


Build and run the test (Linux):
```bash
make
bin/FiDiTi
python plot.py   # Output in './out/plot/'
```

Switch to Release mode (optimized):
```bash
make clean
make TARGET=Release   # Or set TARGET in Makefile
bin/FiDiTi
python plot.py
```



## Installation

To use the library in another project with CMake, install it to a local directory:
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=<Install-dir>
cmake --install build
```

CMake's `find_package()` should then find it when the library location is added to CMake's search paths.



## Integration in ViennaPS

See the `photodiode` example in [ViennaPS](https://github.com/ViennaTools/ViennaPS). It uses this library as an external dependency with CMake.

To enable it, pass the install location of the FDTD library to CMake (variable `VIENNAPS_LOOKUP_DIRS`). It should then find it when configuring the build in ViennaPS and activate the example.
```bash
cmake -B build -DVIENNAPS_BUILD_EXAMPLES=ON -DVIENNAPS_LOOKUP_DIRS=<Install-dir>
cd build/examples/photodiode
make
```
