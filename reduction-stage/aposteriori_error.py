import numpy as np
import sys

targ = np.loadtxt(sys.argv[1])
pred = np.loadtxt(sys.argv[2])

length = np.shape(targ)[0]
length_p = np.shape(pred)[0]

if length_p < length:
        length = length_p

targ = targ[:length, :]
pred = pred[:length, :]

rmse = np.sqrt(((pred - targ)**2).mean(axis=0))
np.set_printoptions(linewidth=120, precision=3)
print(rmse)
