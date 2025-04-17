import numpy as np
from matplotlib import pyplot as plt
import japanize_matplotlib


x = []
eps = 0.1
while(eps > 0.9*1e-25):
    x.append(eps)
    eps = eps / np.sqrt(10)
# 3.16228e-07, 66, 58, 62.66, 2.3842
# 3.16228e-08, 76, 68, 72.8, 2.498
# 3.16228e-10, 96, 88, 92.6, 2.40832
# 3.16228e-09, 86, 76, 82.58, 2.17798
dete_AVE = [5.4800000000000004263, 11.060000000000000497, 16.039999999999999147, 21.019999999999999574, 25.800000000000000711, 30.539999999999999147, 35.74000000000000199, 40.380000000000002558, 45.479999999999996874, 50.520000000000003126, 55.520000000000003126, 60.640000000000000568, 65.780000000000001137, 70.540000000000006253, 75.519999999999996021, 80.439999999999997726, 85.439999999999997726, 90.560000000000002274, 95.620000000000004547, 100.35999999999999943, 105.62000000000000455, 110.29999999999999716, 115.07999999999999829, ]
dete_MIN = [0, 6, 10, 16, 22, 24, 30, 34, 34, 46, 48, 48, 58, 64, 70, 74, 78, 84, 90, 96, 102, 102, 110, ]
dete_MAX = [8, 14, 18, 24, 28, 34, 38, 44, 48, 54, 58, 64, 70, 74, 80, 84, 88, 94, 98, 104, 108, 114, 118, ]
dete_x = x[:len(dete_AVE)]


x2 = []
eps = np.sqrt(0.1)
while(eps >= 1e-10):
    x2.append(eps)
    eps = eps / 10


# MAX2 = [8, 18, 28, 36, 48, 56]
# MIN2 = [0, 8, 18, 26, 38, 48]
# mean2 = [4.46, 13.8, 23.5, 33.28, 43.18, 53.54]
# std2 = [1.71709, 2.40832, 2.38956, 2.63089, 2.38487, 2.16988]
# x2 = x2[:len(MAX2)]


x_3log = [1/(i+1) for i in range(0, int(1e18), int(1e14))]
y_3log = [3*np.log2(1/eps) for eps in x_3log]
y_9log = [9*np.log2(1/eps) for eps in x_3log]
y_15log = [1.5*np.log2(1/eps) for eps in x_3log]


RS_MAX = [48, 64, 82, 94, 108, 128, 142, 156, 170, 190, 206, 220, 236, 250, 264, 278, 292, 310, 326, 336, 352, 366, 382, 396, 412, 426, 442, 456, 472, 488, 504, 518, 532, 546, 562, 580, 592, 608, 624]
RS_MIN = [11, 34, 46, 58, 62, 79, 93, 97, 103, 121, 186, 198, 210, 232, 244, 258, 270, 288, 300, 322, 336, 352, 368, 380, 396, 412, 422, 438, 452, 472, 484, 500, 516, 530, 544, 564, 576, 588, 604]
RS_AVE = [36.63, 52.91, 68.71, 85.23, 101.39, 117.75, 132.87, 147.47, 163.19, 178.23, 194.7, 209.28, 223.94, 239.5, 254.7, 269.66, 283.86, 299.62, 315.18, 330.1, 345.52, 359.82, 375.04, 389.9, 405.1, 420.22, 434.52, 449.8, 464.58, 480.18, 495.24, 510.4, 525.62, 540.36, 555.64, 570.7, 585.24, 600.0, 615.28]
RS_x = x[:len(RS_AVE)]


prob_AVE = [1.7700000000000000178, 3.7599999999999997868, 6.6799999999999997158, 9.1999999999999992895, 11.970000000000000639, 14.050000000000000711, 16.780000000000001137, 18.960000000000000853, 21.609999999999999432, 23.960000000000000853, 26.320000000000000284, 28.769999999999999574, 31.280000000000001137, 33.679999999999999716, 36.320000000000000284, 38.670000000000001705, 41.049999999999997158, 43.820000000000000284, 46.140000000000000568, 48.729999999999996874, 51.469999999999998863, 53.909999999999996589, 56.25, 58.979999999999996874, 61.429999999999999716, 63.789999999999999147, 66.349999999999994316, 68.75, 71.25, 73.870000000000004547, 76.140000000000000568, 78.730000000000003979, 81.180000000000006821, ]
prob_MIN = [0, 2, 4, 5, 9, 11, 14, 16, 19, 21, 23, 24, 29, 31, 33, 34, 37, 41, 43, 45, 48, 48, 53, 56, 58, 60, 63, 65, 67, 69, 71, 76, 78, ]
prob_MAX = [3, 5, 9, 12, 14, 16, 19, 21, 25, 26, 28, 31, 33, 36, 38, 41, 43, 46, 48, 51, 54, 56, 58, 61, 63, 66, 68, 71, 74, 76, 79, 80, 83, ]
prob_x = x[:len(prob_AVE)]

clr = plt.cm.Purples(0.9)
fig, ax = plt.subplots()
ax.scatter(dete_x, dete_AVE, label="決定的 : 平均", color="purple", marker="^")
ax.fill_between(dete_x, dete_MIN, dete_MAX, alpha=0.2, edgecolor="purple", facecolor="purple", label="決定的 : 最小値-最大値", hatch="///")
ax.plot(x_3log, y_3log, "--", color="purple", label="$3$log$_{2}(1/\epsilon)$")
ax.grid()

ax.scatter(RS_x, RS_AVE, label="RS法 : 平均", color="g", marker="^")
ax.fill_between(RS_x, RS_MIN, RS_MAX, alpha=0.2, edgecolor="g", facecolor="g", label="RS法 : 最小値-最大値", hatch="///")
ax.plot(x_3log, y_9log, ":", color="g", label="$9$log$_{2}(1/\epsilon)$")

ax.scatter(prob_x, prob_AVE, label="確率的 : 平均", color="r", marker="^")
ax.fill_between(prob_x, prob_MIN, prob_MAX, alpha=0.2, edgecolor="r", facecolor="r", label="確率的 : 最小値-最大値", hatch="///")
ax.plot(x_3log, y_15log, ":", color="r", label="$1.5$log$_{2}(1/\epsilon)$")

ax.set_xscale("log")
x.reverse()
ax.invert_xaxis()
ax.set_xlim(0.2, 0.5*1e-17)

ax.set_xlabel("近似精度 ε", fontsize=18)
ax.set_ylabel("Tゲート数", fontsize=18)
# ax.set_title("ランダムユニタリの$\epsilon$近似に必要なTゲート数", fontsize=18)

ax.set_xticks([10**(-i) for i in range(1, 18, 2)])

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

plt.savefig("plot/T_count.pdf", dpi=500)
