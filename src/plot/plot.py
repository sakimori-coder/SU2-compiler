import numpy as np
from matplotlib import pyplot as plt


x = []
eps = 0.1
while(eps >= 1e-10):
    x.append(eps)
    eps = eps / 10

MAX = [12, 22, 30, 40, 52, 60, 70]
MIN = [2, 10, 20, 30, 46, 54, 58]
mean = [7.96, 18, 27.78, 38.04, 48.04, 57.92, 67.7]
x = x[:len(MAX)]

x_ideal = [1/(i+1) for i in range(0, int(1e11), 100000)]
y_ideal = [3*np.log2(1/eps) for eps in x_ideal]

# x.reverse()
# MAX.reverse()
# MIN.reverse()
# mean.reverse()

clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(x, mean)
ax.fill_between(x, MIN, MAX, alpha=0.3, edgecolor=clr, facecolor=clr)
ax.plot(x_ideal, y_ideal)

ax.set_xscale("log")
x.reverse()
ax.set_xticks(x)
ax.invert_xaxis()
ax.set_xlim(0.5, 1e-8)

plt.savefig("T_count.png")
