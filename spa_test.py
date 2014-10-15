# -*- coding: utf-8 -*-
import numpy as np
import _spa

yr=np.array([[1988, 1988], [2002, 2003], [2010, 2014]])
mo=np.array([[6, 12], [5, 4], [1, 12]])
da=np.array([[1, 20], [22, 3], [5, 5]])
ho=np.array([[0, 23], [10, 5], [0, 22]])
mi=np.array([[23, 59], [0, 33], [1, 1]])
se=np.array([[23, 59], [0, 3], [3, 3]])
la=np.array([[10, -20], [-90, 90], [-90, 90]])
lo=np.array([[100, 270], [-130, 300], [-180, 180]])
el=np.array([[0, -100], [30, 300], [3000, 5000]])
sl=np.array([[0, -10], [1, 40], [0, 6]])
ap=np.array([[100, 270], [-130, -300], [0, 330]])

ze,az=_spa.za(yr, mo, da, ho, mi, se, la, lo)

print ze
print az

ze,az,ic=_spa.za_inc(yr, mo, da, ho, mi, se, la, lo, el, sl, ap)

print 'za_inc'
print ze
print ic

# pygrads test

import grads, os
import matplotlib.pyplot as plt

gradsbase="/home/"+os.environ["localbase"]+"/"+os.environ["USER"]+"/local/opengrads"
ga = grads.GrADS(Bin=gradsbase+"/Linux/x86_64/grads", Window=False)

ga.open(gradsbase+"/Resources/SampleDatasets/model.ctl")

dummy=ga.exp('ts')
zz=dummy.data*0

yr=zz+2003
mo=zz+12
da=zz+31
ho=zz+23
mi=zz+59
se=zz+59
lo,la=np.meshgrid(dummy.grid.lon, dummy.grid.lat)

ze,az=_spa.za(yr, mo, da, ho, mi, se, la, lo)

dummy.data[:]=ze[:]
ga.imp('ze', dummy)

ga.basemap('cyl')
ga.contourf('ze')
plt.title('Solar Zenith Angle')
plt.savefig('test_za.png')
