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


x_3log = [1/(i+1) for i in range(0, int(1e16), int(1e14))]
y_3log = [3*np.log2(1/eps) for eps in x_3log]
y_9log = [9*np.log2(1/eps) for eps in x_3log]
y_15log = [1.5*np.log2(1/eps) for eps in x_3log]

Rs_x = [np.power(np.sqrt(0.1), i) for i in range(1, 21)]
RS_MAX = [3, 26, 48, 64, 78, 94, 110, 126, 142, 158, 172, 186, 202, 216, 232, 246, 264, 276, 296, 312]
RS_MIN = [0, 0, 12, 21, 49, 53, 90, 83, 124, 142, 154, 162, 184, 200, 216, 230, 248, 262, 274, 292]
RS_mean = [1.34, 12.07, 36.05, 54.53, 70.85, 85.53, 103.2, 118.23, 133.9, 149.48, 164.72, 179.22, 194.18, 209.4, 225.38, 239.8, 254.92, 269.46, 285.48, 300.76]

prob_eps = [0.10000000000000000555, 0.031622776601683795873, 0.01000000000000000106, 0.0031622776601683797468, 0.0010000000000000001564, 0.00031622776601683799064, 0.00010000000000000002069, 3.1622776601683800659e-05, 1.0000000000000002573e-05, 3.1622776601683802254e-06, 1.0000000000000003078e-06, 3.162277660168380385e-07, 1.0000000000000003582e-07, 3.1622776601683805445e-08, 1.0000000000000004086e-08, 3.162277660168380704e-09, 1.0000000000000004591e-09, 3.1622776601683808636e-10, 1.0000000000000005095e-10] + [3.162277660168379527e-11, 1.0000000000000000869e-11, 3.1622776601683796865e-12, 1.0000000000000001373e-12, 3.162277660168379846e-13, 1.0000000000000001878e-13, 3.1622776601683800056e-14, 1.0000000000000002382e-14, 3.1622776601683801651e-15, 1.0000000000000002887e-15]
prob_AVE = [1.5300000000000000266, 3.4599999999999999645, 7.1500000000000003553, 10, 13.099999999999999645, 15.019999999999999574, 17.460000000000000853, 20.039999999999999147, 22.640000000000000568, 24.920000000000001705, 27.339999999999999858, 30.239999999999998437, 32.579999999999998295, 34.920000000000001705, 37.359999999999999432, 39.799999999999997158, 42.600000000000001421, 45.039999999999999147, 47.520000000000003126] + [49.979999999999996874, 52.619999999999997442, 54.960000000000000853, 57.399999999999998579, 59.74000000000000199, 62.359999999999999432, 65, 67.280000000000001137, 69.680000000000006821, 72.560000000000002274, ]
prob_MIN = [0, 1, 4, 3, 10, 12, 14, 16, 20, 20, 24, 28, 30, 32, 34, 36, 40, 42, 44,] + [48, 48, 50, 54, 52, 60, 62, 64, 66, 70, ]
prob_MAX = [2, 5, 10, 12, 16, 18, 20, 22, 24, 28, 30, 32, 34, 38, 40, 42, 44, 48, 50] + [52, 56, 58, 60, 62, 64, 66, 70, 72, 74, ]


clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(x, mean, label="Direct : ave", color="purple")
ax.fill_between(x, MIN, MAX, alpha=0.2, edgecolor="purple", facecolor="purple", label="Direct : min-max")
ax.plot(x_3log, y_3log, "--", color="purple", label="$3$log$_{2}(1/\epsilon)$")
ax.grid()

ax.scatter(Rs_x, RS_mean, label="RS : ave", color="g", marker="^")
ax.fill_between(Rs_x, RS_MIN, RS_MAX, alpha=0.2, edgecolor="g", facecolor="g", label="RS : min-max", hatch="///")
ax.plot(x_3log, y_9log, ":", color="g", label="$9$log$_{2}(1/\epsilon)$")

ax.scatter(prob_eps, prob_AVE, label="Mixed : ave", color="r", marker="^")
ax.fill_between(prob_eps, prob_MIN, prob_MAX, alpha=0.2, edgecolor="r", facecolor="r", label="Mixed : min-max", hatch="///")
ax.plot(x_3log, y_15log, ":", color="r", label="$1.5$log$_{2}(1/\epsilon)$")

ax.set_xscale("log")
x.reverse()
ax.set_xticks(prob_eps)
ax.invert_xaxis()
ax.set_xlim(0.2, 0.5*1e-15)

ax.set_xlabel("Approximation accuracy ε", fontsize=18)
ax.set_ylabel("T-count", fontsize=18)
# ax.set_title("ランダムユニタリの$\epsilon$近似に必要なTゲート数", fontsize=18)

handles, labels = ax.get_legend_handles_labels()
new_handles = []
new_labels = []
new_handles.append(handles[3])
new_handles.append(handles[4])
new_handles.append(handles[0])

new_handles.append(handles[7])
new_handles.append(handles[8])
new_handles.append(handles[2])

new_handles.append(handles[5])
new_handles.append(handles[6])
new_handles.append(handles[1])



new_labels.append(labels[3])
new_labels.append(labels[4])
new_labels.append(labels[0])

new_labels.append(labels[7])
new_labels.append(labels[8])
new_labels.append(labels[2])

new_labels.append(labels[5])
new_labels.append(labels[6])
new_labels.append(labels[1])

ax.legend(loc="upper left", handles=new_handles, labels=new_labels)

plt.savefig("T_count.png", dpi=500)
