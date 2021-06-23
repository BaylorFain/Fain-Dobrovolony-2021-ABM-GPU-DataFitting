#README: This code compares the run times of a code written in three different ways

import numpy
import datetime
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import anderson
from scipy.stats import normaltest
from scipy.stats import ttest_1samp
import scipy.integrate as integrate
import matplotlib.ticker as ticker

C2      = [0.175063, 0.178323, 0.179585, 0.179549, 0.174793, 0.178617, 0.175439, 0.172861, 0.168731, 0.163323]
C3      = [1.4927, 1.46637, 1.46689, 1.47239, 1.48267, 1.46128, 1.46301, 1.46176, 1.478, 1.46468]
C4      = [15.3137, 15.2461, 15.1718, 15.1554, 15.2294, 15.2685, 15.1771, 15.2608, 15.2902, 15.243]
C5      = [159.345, 154.286, 154.156, 154.237, 154.452, 154.183, 154.414, 154.043, 154.711, 153.594]
C6      = [535.781, 489.157, 495.174, 489.615, 489.172, 488.846, 487.327, 488.204, 487.875, 488.401]

CUDA2   = [0.0292391, 0.0256742, 0.0249308, 0.0251271, 0.025154, 0.0252087, 0.0249681, 0.0254071, 0.0255272, 0.023015]
CUDA3   = [0.0335051, 0.0343513, 0.0334941, 0.0335641, 0.0334063, 0.0313023, 0.0308947, 0.0310007, 0.0311123, 0.0313072]
CUDA4   = [0.0766135, 0.0740529, 0.0753639, 0.0704517, 0.0704176, 0.0651617, 0.0652225, 0.0660366, 0.0618907, 0.0633054]
CUDA5   = [0.565216, 0.567913, 0.543877, 0.515357, 0.51727, 0.525859, 0.518283, 0.51908, 0.525726, 0.514576]
CUDA6   = [11.377, 11.4964, 11.494, 11.4581, 11.4539, 11.414, 11.6031, 11.4949, 11.5166, 11.6028]

Python2 = [12.199251, 11.633082, 11.656383, 11.568248, 11.680383, 11.667641, 11.622523, 11.605636, 11.689679, 11.649330]
Python3 = [54.454668, 55.487442, 55.504977, 55.550025, 55.309193, 55.844685, 55.226034, 55.458564, 55.356308, 55.315955]
Python4 = [477.393221, 488.609290, 489.176073, 487.582017, 489.742517, 490.387996, 490.467258, 489.883759, 489.642600, 489.744665]
Python5 = [5334.968568, 5468.218875, 5463.220027, 5488.610683, 5482.253036, 5498.881454, 5498.278449, 5504.448134, 5500.209880, 5491.251355]
Python6 = [78089.182268, 80609.944895, 80597.736797, 80902.591233, 80708.631306, 79819.988632, 79549.839675, 79848.759659, 81278.175224, 77610.220007]

C = []
CUDA = []
Python = []

C.append(numpy.mean(C2))
C.append(numpy.mean(C3))
C.append(numpy.mean(C4))
C.append(numpy.mean(C5))
C.append(numpy.mean(C6))

CUDA.append(numpy.mean(CUDA2))
CUDA.append(numpy.mean(CUDA3))
CUDA.append(numpy.mean(CUDA4))
CUDA.append(numpy.mean(CUDA5))
CUDA.append(numpy.mean(CUDA6))

Python.append(numpy.mean(Python2))
Python.append(numpy.mean(Python3))
Python.append(numpy.mean(Python4))
Python.append(numpy.mean(Python5))
Python.append(numpy.mean(Python6))

print(C)
print(CUDA)
print(Python)

NumberofCells = [121, 1027, 10009, 100969, 1001365]

size = 32
plt.rcParams['xtick.labelsize'] = size
plt.rcParams['ytick.labelsize'] = size
plt.rcParams['axes.labelsize'] = size
plt.rcParams['figure.figsize'] = [11, 11]
fig, ax = plt.subplots()
ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

plt.grid(b="TRUE", axis="both")
plt.xscale("log")
plt.yscale("log")
plt.plot(NumberofCells, numpy.divide(C,1),"o-", linewidth = 3, label = "C")
plt.plot(NumberofCells, numpy.divide(CUDA,1),"o-", linewidth = 3, label = "CUDA")
plt.plot(NumberofCells, numpy.divide(Python,1),"o-", linewidth = 3, label = "Python")
plt.xlabel("Number of cells")
plt.ylabel("Computation time of one simulated hour (s)")
plt.legend(prop={'size': size*(7/10)})
plt.tight_layout()
plt.close()
fig.savefig(os.path.join("/home/baylor/Documents/Research/Papers/Current_Papers/Computer_Science/Figures","loglogspeed"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
