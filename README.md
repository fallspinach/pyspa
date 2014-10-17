pyspa
=====

Python module for the Solar Position Algorithm

This is a python wrapper module for the Solar Position Algorithm (SPA) written by Afshin Michael Andreas <Afshin.Andreas@NREL.gov>.

Ming Pan <mpan@Princeton.EDU>

###Build

python setup.py build_ext --inplace

###Usage

The module provides one function _spa.calc().

```python
import _spa

zenith,azimuth = _spa.calc(year, month, day, hour, minute, second, latitude, longitude)

zenith,azimuth,incidence = _spa.calc(year, month, day, hour, minute, second, latitude, longitude, elevation, slope, aspect)

```

The function requires either exactly 8 for zenith/azimuth calculation or exactly 11 inputs for zenith/azimuth/incidence calculation. All inputs can either be a single value varialbe or a numpy array. There must be at least one numpy array (even with just one element) and all the arrays must have the same dimensions.

###Test

python spa_test.py

**Note**: pygrads with opengrads is needed to run the second part of the test.

