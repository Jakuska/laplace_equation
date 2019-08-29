import easygui
from matplotlib import pyplot as plt
from matplotlib import image as image
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from video_function import *

N = 60
w, h = N/3, (2*N)/3
prism = np.zeros((N, N), dtype=float)
boundary_voltage = 1
boundary_location = []
for i in range(w, h):
    for j in range(w, h):
            prism[i][j] = boundary_voltage
            boundary_location.append((i,j))

j_result = two_dimension_jacobi(prism, boundary_location, 0.01)

print(len(j_result[3]))

for i in range(len(j_result[3])):
    matrix = np.array(j_result[3][i])
    xlist = np.linspace(-1.0, 1.0, N)
    ylist = np.linspace(-1.0, 1.0, N)
    X, Y = np.meshgrid(xlist, ylist)
    Z = matrix[0]

    fig = plt.figure(figsize=(20, 20))
    plt.axis('off')
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0,
                    cmap='bone', edgecolor='none');

    plt.savefig('video_images/prism_%d.png' % i)
    print(i)
    fig.clf()
    plt.close('all')
