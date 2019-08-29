import easygui
from matplotlib import pyplot as plt
from matplotlib import image as image
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numerical_methods import *

width, height = 30, 30
offset = 5
matrix_midpoint = width/2
capacitor = np.zeros((width, height), dtype=float)
boundary_location = []
for row in range(int(height * 0.3), int(height * 0.7)):
    boundary_location.append((row, matrix_midpoint - offset))
    boundary_location.append((row, matrix_midpoint + offset))
    capacitor[row][matrix_midpoint - offset] = -1
    capacitor[row][matrix_midpoint + offset] = +1


output_matrix = two_dimension_jacobi(capacitor, boundary_location, microVolt, middle_matrix=25)

fig, (ax1,ax2,ax3) = plt.subplots(1,3)
first = ax1.imshow(capacitor, cmap="bone")
first = ax2.imshow(output_matrix[1], cmap="bone")
first = ax3.imshow(output_matrix[0], cmap="bone")
fig.colorbar(first, ax=[ax1, ax2, ax3], orientation="horizontal")

plt.show()

# Show Voltage
# jacobi_solution = output_matrix[0]
# coordinates  = np.linspace(0, 1.0, len(jacobi_solution))
# X, Y   = np.meshgrid(coordinates, coordinates)
# dx, dy = np.diff(X, axis=1).item(0), np.diff(Y, axis=0).item(0)
# Ey, Ex = np.gradient(- jacobi_solution)
# Z      = np.sqrt((Ex/dx) ** 2 + (Ey/dy) **2)
# np.savetxt("Z_capacitor.csv", Z, delimiter=",")
# plt.streamplot(X, Y, Ex, Ey, density=1, color='k', linewidth=Z/5)

# plt.show()
