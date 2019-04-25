import numpy
import matplotlib.pyplot as plt

d_compression = numpy.loadtxt("strain-stress-damage_compression.dat")[:, [0, 6]]
d_traction = numpy.loadtxt("strain-stress-damage_traction.dat")[:, [0, 6]]
dto_compression = numpy.loadtxt("strain-stress-damage-traction-only_compression.dat")[:, [0, 6]]
dto_traction = numpy.loadtxt("strain-stress-damage-traction-only_traction.dat")[:, [0, 6]]

plt.figure(1)

plt.subplot(221)
plt.title("Damage - traction")
plt.plot(d_traction[:,0], d_traction[:, 1])
plt.ylabel("Stress XX")

plt.subplot(222)
plt.title("Damage - compression")
plt.plot(d_compression[:,0], d_compression[:, 1])

plt.subplot(223)
plt.title("DamageTractionOnly - traction")
plt.plot(dto_traction[:,0], dto_traction[:, 1])
plt.ylabel("Stress XX")
plt.xlabel("Strain XX")

plt.subplot(224)
plt.title("DamageTractionOnly - compression")
plt.plot(dto_compression[:,0], dto_compression[:, 1])
plt.xlabel("Strain XX")

plt.show()
