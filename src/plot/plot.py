import numpy as np
from matplotlib import pyplot as plt
import japanize_matplotlib


x = []
eps = np.sqrt(0.1)
while(eps > 0.9*1e-10):
    x.append(eps)
    eps = eps / np.sqrt(10)
# 3.16228e-07, 66, 58, 62.66, 2.3842
# 3.16228e-08, 76, 68, 72.8, 2.498
# 3.16228e-10, 96, 88, 92.6, 2.40832
# 3.16228e-09, 86, 76, 82.58, 2.17798
MAX = [8, 12, 18, 22, 28, 30, 36, 40, 48, 52, 56, 60, 66, 70, 76, 80, 86, 90, 96, 100]
MIN = [0, 2, 8, 10, 18, 20, 26, 30, 38, 46, 48, 54, 58, 58, 68, 68, 76, 78, 88, 90]
mean = [4.46, 7.96, 13.8, 18, 23.5, 27.78, 33.28, 38.04, 43.18, 48.04, 53.54, 57.92, 62.66, 67.7, 72.8, 77.52, 82.58, 87.86, 92.6, 98]
x = x[:len(MAX)]
print(x)
print(len(MAX))

x2 = []
eps = np.sqrt(0.1)
while(eps >= 1e-10):
    x2.append(eps)
    eps = eps / 10


MAX2 = [8, 18, 28, 36, 48, 56]
MIN2 = [0, 8, 18, 26, 38, 48]
mean2 = [4.46, 13.8, 23.5, 33.28, 43.18, 53.54]
std2 = [1.71709, 2.40832, 2.38956, 2.63089, 2.38487, 2.16988]
x2 = x2[:len(MAX2)]


x_3log = [1/(i+1) for i in range(0, int(1e11), 10000000)]
y_3log = [3*np.log2(1/eps) for eps in x_3log]
y_9log = [9*np.log2(1/eps) for eps in x_3log]


Rs_x = [np.power(np.sqrt(0.1), i) for i in range(1, 21)]
RS_MAX = [3, 26, 48, 64, 78, 94, 110, 126, 142, 158, 172, 186, 202, 216, 232, 246, 264, 276, 296, 312]
RS_MIN = [0, 0, 12, 21, 49, 53, 90, 83, 124, 142, 154, 162, 184, 200, 216, 230, 248, 262, 274, 292]
RS_mean = [1.34, 12.07, 36.05, 54.53, 70.85, 85.53, 103.2, 118.23, 133.9, 149.48, 164.72, 179.22, 194.18, 209.4, 225.38, 239.8, 254.92, 269.46, 285.48, 300.76]


clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(x, mean, label="提案 : 平均値", color="purple")
ax.fill_between(x, MIN, MAX, alpha=0.3, edgecolor="purple", facecolor="purple", label="提案 : 最小値-最大値")
ax.plot(x_3log, y_3log, "--", color="purple", label="$3log_{2}(1/\epsilon)$")
ax.grid()

ax.scatter(Rs_x, RS_mean, label="RS : 平均値", color="g")
ax.fill_between(Rs_x, RS_MIN, RS_MAX, alpha=0.3, edgecolor="g", facecolor="g", label="RS : 最小値-最大値")
ax.plot(x_3log, y_9log, "--", color="g", label="$9log_{2}(1/\epsilon)$")


ax.set_xscale("log")
x.reverse()
ax.set_xticks(x)
ax.invert_xaxis()
ax.set_xlim(0.7, 0.5*1e-10)

ax.set_xlabel("Approximation accuracy ε", fontsize=18)
ax.set_ylabel("T-count", fontsize=18)
# ax.set_title("ランダムユニタリの$\epsilon$近似に必要なTゲート数", fontsize=18)

handles, labels = ax.get_legend_handles_labels()
new_handles = []
new_labels = []
new_handles.append(handles[2])
new_handles.append(handles[3])
new_handles.append(handles[0])
new_handles.append(handles[4])
new_handles.append(handles[5])
new_handles.append(handles[1])

new_labels.append(labels[2])
new_labels.append(labels[3])
new_labels.append(labels[0])
new_labels.append(labels[4])
new_labels.append(labels[5])
new_labels.append(labels[1])

ax.legend(loc="upper left", handles=new_handles, labels=new_labels)

plt.savefig("T_count.png", dpi=500)
