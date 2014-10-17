pyspa
=====

Python module for the Solar Position Algorithm

This is a python wrapper module for the Solar Position Algorithm (SPA) written by Afshin Michael Andreas <Afshin.Andreas@NREL.gov>.

Ming Pan <mpan@Princeton.EDU>

###Build

python setup.py build_ext --inplace

###Usage

import _spa

ze,az = _spa.calc(yr, mo, da, ho, mi, se, la, lo)

ze,az,ic = _spa.calc(yr, mo, da, ho, mi, se, la, lo, el, sl, ap)

The function requires either exactly 8 for zenith/azimuth calculation or exactly 11 inputs for zenith/azimuth/incidence calculation. All inputs can either be a single value varialbe or a numpy array. There must be at least one numpy array (even with just one element) and all the arrays must have the same dimensions.

###Test

python spa_test.py

**Note**: pygrads with opengrads is needed to run the second part of the test.

