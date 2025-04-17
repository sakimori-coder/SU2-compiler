import numpy as np
from matplotlib import pyplot as plt
import japanize_matplotlib



eps = [0.10000000000000000555, 0.031622776601683795873, 0.01000000000000000106, 0.0031622776601683797468, 0.0010000000000000001564, 0.00031622776601683799064, 0.00010000000000000002069, 3.1622776601683800659e-05, 1.0000000000000002573e-05, 3.1622776601683802254e-06, 1.0000000000000003078e-06, 3.162277660168380385e-07, 1.0000000000000003582e-07, 3.1622776601683805445e-08, 1.0000000000000004086e-08, 3.162277660168380704e-09, 1.0000000000000004591e-09, 3.1622776601683808636e-10, 1.0000000000000005095e-10, 3.1622776601683810231e-11, 1.00000000000000056e-11, 3.1622776601683811826e-12, 1.0000000000000006104e-12, ]
TIME = [6.4382700000000019358, 17.825220000000005172, 34.780830000000001689, 55.327199999999983504, 81.187819999999987886, 113.25453999999994892, 154.95674999999991428, 218.83597999999997796, 281.4052499999999668, 364.18639000000030137, 475.96327999999999747, 639.84685999999976502, 919.06865999999979522, 1326.6061000000001968, 1974.8718799999999192, 3091.8330399999990732, 5450.51987999999983, 10360.130230000002484, 25253.644619999988208, 57925.088300000003073, 128604.283310000028, 287948.18850000010571, 655588.96992000017781, ]
eps = eps[:len(TIME)]


print(len(eps))
print(len(TIME))
print(eps[:len(TIME)+1])

log_eps = np.log(eps)
log_TIME = np.log(TIME)
a, b = np.polyfit(log_eps[-10:], log_TIME[-10:], 1)
print(a, np.exp(b))

# y_scale = [0.3 * np.power(eps[i], -0.6) for i in range(len(eps))]
y_scale = [np.exp(b) * np.power(eps[i], a) for i in range(len(eps))]


clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(eps, TIME, color="red")
ax.plot(eps, y_scale, "--", color="blue", label="$t = 0.019\epsilon^{-0.62}$")
ax.grid()


ax.set_xscale("log")
ax.set_yscale("log")
eps.reverse()
ax.set_xticks(eps)
ax.invert_xaxis()
# ax.set_xlim(0.2, 0.5*1e-10)
ax.set_ylim(1e0, 1e6)

ax.set_xlabel("近似精度 ε", fontsize=18)
ax.set_ylabel("時間[ms]", fontsize=18)
# ax.set_title("決定的Clifford+$T$分解アルゴリズム", fontsize=18)

ax.set_xticks([10**(-i) for i in range(1, 13, 1)])

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
ax.legend(fontsize=15)

plt.savefig("plot/dete_time.pdf", dpi=500)
