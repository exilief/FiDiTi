
<h1>FiDiTi - Finite-Difference Time-Domain Method (FDTD)</h1>

Header-only C++ library


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
