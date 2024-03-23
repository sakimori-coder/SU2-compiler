from matplotlib import pyplot as plt
# import japanize_matplotlib
import matplotlib.patches as patches
import matplotlib.animation as animation
import random
import math
import copy


def distance(xy1, xy2):
    return math.sqrt((xy1[0] - xy2[0])**2 + (xy1[1] - xy2[1])**2)

def cal_inter(xy1, xy2, r1, r2):
    x1 = xy1[0]
    y1 = xy1[1]
    x2 = xy2[0]
    y2 = xy2[1]
    a = 2 * (x2 - x1)
    b = 2 * (y2 - y1)
    c = (x1 + x2)*(x1 - x2) + (y1 + y2)*(y1 - y2) + (r2 + r1)*(r2 - r1) 

    D = -(a*x1 + b*y1 + c)
    E = (a**2 + b**2) * (r1**2) - D**2

    if E < 0:
        return None
    
    Ax = (a*D - b*math.sqrt(E)) / (a**2 + b**2) + x1
    Ay = (b*D + a*math.sqrt(E)) / (a**2 + b**2) + y1
    
    Bx = (a*D + b*math.sqrt(E)) / (a**2 + b**2) + x1
    By = (b*D - a*math.sqrt(E)) / (a**2 + b**2) + y1

    return (Ax, Ay), (Bx, By)



fig = plt.figure()
ax = plt.axes()

marker1 = '😁'   # 成功
marker2 = '😰'   # 失敗
marker3 = '😐'   # 考え中
emoji1 = ax.annotate(marker1, (2.5, 2.5),fontsize=50, ha='center',va='center',
            xytext=(0,0), textcoords='offset pixels')
emoji2 = ax.annotate(marker2, (2.5, 2.5),fontsize=50, ha='center',va='center',
            xytext=(0,0), textcoords='offset pixels')
emoji3 = ax.annotate(marker3, (2.5, 2.5),fontsize=50, ha='center',va='center',
            xytext=(0,0), textcoords='offset pixels')

eps = 1
r_large = 2*eps
large_ball = patches.Circle(xy=(0, 0), radius=r_large, fc='g')
ax.add_patch(large_ball)

random.seed(5)

n = 20   # 点の数

xy = []
r_small = eps
for i in range(n):
    x = random.uniform(-r_large, r_large)
    y = random.uniform(-math.sqrt(r_large**2 - x**2), math.sqrt(r_large**2 - x**2))
    xy.append((x,y))

ball_list = []
for j in range(n):
    x = xy[j][0]
    y = xy[j][1]
    ball = patches.Circle(xy=(x,y), radius=r_small, ec='b', fill=False)
    ax.add_patch(ball)
    ball_list.append(ball)

red_balls = []
for j in range(n):
    x = xy[j][0]
    y = xy[j][1]
    ball = patches.Circle(xy=(x,y), radius=r_small, ec='r', fill=False)
    # ax.add_patch(ball)
    red_balls.append(ball)

animations = []

for i in range(n):
    for j in range(i+1,n):
        ret = cal_inter(xy[i], xy[j], r_small, r_small)
        if ret == None:
            continue
        A, B = ret

        if distance(A, (0,0)) < r_large:
            s1 = ax.scatter(A[0], A[1], color='r')
            animations.append([s1] + ball_list + [emoji3])
            flag = False
            for k in range(n):
                if k == i or k == j:
                    continue
                if distance(A, xy[k]) < r_small:
                    art = ax.add_artist(red_balls[k])
                    animations.append([art] + ball_list + [s1] + [emoji1])
                    flag = True
                    break
            if not flag: 
                animations.append([s1] + ball_list + [emoji2])
                # exit(0)

        if distance(B, (0,0)) < r_large:
            s2 = ax.scatter(B[0], B[1], color='r')
            animations.append([s2] + ball_list + [emoji3])
            flag = False
            for k in range(n):
                if k == i or k == j:
                    continue
                if distance(B, xy[k]) < r_small:
                    art = ax.add_artist(red_balls[k])
                    animations.append([art] + ball_list + [s2] + [emoji1])
                    flag = True
                    break
            if not flag: 
                animations.append([s2] + ball_list + [emoji2])



for i in range(n):
    ret = cal_inter(xy[i], (0,0), r_small, r_large)
    if ret == None:
        continue
    A, B = ret

    s1 = ax.scatter(A[0], A[1], color='r')
    animations.append([s1] + ball_list + [emoji3])
    flag = False
    for j in range(n):
        if j == i:
            continue
        if distance(A, xy[j]) < r_small:
            art = ax.add_artist(red_balls[j])
            animations.append([art] + ball_list + [s1] + [emoji1])
            flag = True
            break
    if not flag:
        animations.append([s1] + ball_list + [emoji2])

    s2 = ax.scatter(B[0], B[1], color='r')
    animations.append([s2] + ball_list + [emoji3])
    flag = False
    for j in range(n):
        if j == i:
            continue
        if distance(B, xy[j]) < r_small:
            art = ax.add_artist(red_balls[j])
            animations.append([art] + ball_list + [s2] + [emoji1])
            flag = True
            break
    if not flag:
        animations.append([s2] + ball_list + [emoji2])

    


plt.axis('scaled')
ax.set_aspect('equal')

ani = animation.ArtistAnimation(fig, animations, interval=500)
ani.save('algo.gif', writer="imagemagick")
