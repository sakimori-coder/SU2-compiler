import numpy as np
from matplotlib import pyplot as plt
import japanize_matplotlib


x = []
eps = 0.1
while(eps >= 1e-10):
    x.append(eps)
    eps = eps / 10

MAX = [12, 22, 30, 40, 52, 60, 70, 80, 90]
MIN = [2, 10, 20, 30, 46, 54, 58, 68, 78]
mean = [7.96, 18, 27.78, 38.04, 48.04, 57.92, 67.7, 77.52, 87.86]
x = x[:len(MAX)]

x_ideal = [1/(i+1) for i in range(0, int(1e11), 10000000)]
y_ideal = [3*np.log2(1/eps) for eps in x_ideal]

# x.reverse()
# MAX.reverse()
# MIN.reverse()
# mean.reverse()

clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(x, mean, label="平均値", color="b")
ax.fill_between(x, MIN, MAX, alpha=0.3, edgecolor=clr, facecolor=clr, label="最小値-最大値")
ax.plot(x_ideal, y_ideal, "--", color="r", label="$3log_{2}(1/\epsilon)$")
ax.grid()

ax.set_xscale("log")
x.reverse()
ax.set_xticks(x)
ax.invert_xaxis()
ax.set_xlim(0.7, 0.5*1e-9)

ax.set_xlabel("Approximation accuracy ε", fontsize=18)
ax.set_ylabel("T-count", fontsize=18)
ax.set_title("ランダムユニタリの$\epsilon$近似に必要なTゲート数", fontsize=18)

ax.legend(loc="upper left", fontsize=15)

plt.savefig("T_count.png", dpi=500)
