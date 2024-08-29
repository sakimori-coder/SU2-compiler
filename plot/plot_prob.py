import numpy as np
from matplotlib import pyplot as plt
import japanize_matplotlib



eps = [0.10000000000000000555, 0.031622776601683795873, 0.01000000000000000106, 0.0031622776601683797468, 0.0010000000000000001564, 0.00031622776601683799064, 0.00010000000000000002069, 3.1622776601683800659e-05, 1.0000000000000002573e-05, 3.1622776601683802254e-06, 1.0000000000000003078e-06, 3.162277660168380385e-07, 1.0000000000000003582e-07, 3.1622776601683805445e-08, 1.0000000000000004086e-08, 3.162277660168380704e-09, 1.0000000000000004591e-09, 3.1622776601683808636e-10, 1.0000000000000005095e-10, ]
AVE = [1.5300000000000000266, 3.4599999999999999645, 7.1500000000000003553, 10, 13.099999999999999645, 15.019999999999999574, 17.460000000000000853, 20.039999999999999147, 22.640000000000000568, 24.920000000000001705, 27.339999999999999858, 30.239999999999998437, 32.579999999999998295, 34.920000000000001705, 37.359999999999999432, 39.799999999999997158, 42.600000000000001421, 45.039999999999999147, 47.520000000000003126, ]
MIN = [0, 1, 4, 3, 10, 12, 14, 16, 20, 20, 24, 28, 30, 32, 34, 36, 40, 42, 44, ]
MAX = [2, 5, 10, 12, 16, 18, 20, 22, 24, 28, 30, 32, 34, 38, 40, 42, 44, 48, 50, ]

x_15log = [1/(i+1) for i in range(0, int(1e11), 10000000)]
y_15log = [1.5*np.log2(1/eps) for eps in x_15log]


# RS_MAX = [26, 48, 62, 78, 98, 112, 128, 140, 154, 174, 186, 200, 220, 232, 248, 262, 276, 294, 304, 322]
# RS_MIN = [0, 10, 15, 46, 72, 94, 104, 120, 138, 154, 172, 180, 198, 216, 230, 244, 256, 272, 288, 304]
# RS_mean = [11.78, 33.84, 53.44, 69.65, 86.1, 102.54, 117.64, 132.34, 148.38, 163.34, 178.42, 193.48, 208.48, 224.66, 239.06, 254.04, 269.0, 284.86, 299.42, 314.42]


clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(eps, AVE, label="Propose : ave", color="purple")
ax.fill_between(eps, MIN, MAX, alpha=0.2, edgecolor="purple", facecolor="purple", label="Propose : min-max")
ax.plot(x_15log, y_15log, "--", color="purple", label="$3$log$_{2}(1/\epsilon)$")
ax.grid()


ax.set_xscale("log")
eps.reverse()
ax.set_xticks(eps)
ax.invert_xaxis()
ax.set_xlim(0.2, 0.5*1e-10)

ax.set_xlabel("Approximation accuracy ε", fontsize=18)
ax.set_ylabel("T-count", fontsize=18)
# ax.set_title("ランダムユニタリの$\epsilon$近似に必要なTゲート数", fontsize=18)

# handles, labels = ax.get_legend_handles_labels()
# new_handles = []
# new_labels = []
# new_handles.append(handles[2])
# new_handles.append(handles[3])
# new_handles.append(handles[0])
# new_handles.append(handles[4])
# new_handles.append(handles[5])
# new_handles.append(handles[1])

# new_labels.append(labels[2])
# new_labels.append(labels[3])
# new_labels.append(labels[0])
# new_labels.append(labels[4])
# new_labels.append(labels[5])
# new_labels.append(labels[1])

# ax.legend(loc="upper left", handles=new_handles, labels=new_labels)

plt.savefig("plot/prob_T_count.png", dpi=500)
