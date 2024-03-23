# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# fig = plt.figure()

# ims = []

# for i in range(10):
#         rand = np.random.randn(100)     # 100個の乱数を生成
#         im = plt.scatter(rand[0], rand[1])
#         # im = plt.plot(rand)             # 乱数をグラフにする
#         ims.append([im])                  # グラフを配列 ims に追加

# # 10枚のプロットを 100ms ごとに表示
# ani = animation.ArtistAnimation(fig, ims, interval=100)
# plt.show()



import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

nx = 20
ny = 20

fig = plt.figure()
plt.axis([0,nx,0,ny])
ax = plt.gca()
ax.set_aspect(1)

def init():
    # initialize an empty list of cirlces
    return []

def animate(i):
    # draw circles, select to color for the circles based on the input argument i. 
    someColors = ['r', 'b', 'g', 'm', 'y']
    patches = []
    for x in range(0,nx):
        for y in range(0,ny):
            patches.append(ax.add_patch( plt.Circle((x+0.5,y+0.5),0.45,color=someColors[i % 5]) ))
    return patches

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=10, interval=20, blit=True)
plt.show()