from matplotlib import pyplot as plt
from matplotlib import image as image
import numpy as np
from numerical_methods import *

width, height = 20, 20
minima, maxima = -1.0, 1.0 
metal_box = np.zeros((width, height), dtype=float)

boundary_location = []
for i in range(width):                      
    metal_box[i][0] = minima 
    metal_box[i][width-1] = maxima
    metal_box[0][i] = minima + ((maxima - minima) / (width-1) * i)
    metal_box[width-1][i] = minima + ((maxima - minima) / (width-1) * i)
    boundary_location.append((i,0))
    boundary_location.append((i,width-1))
    boundary_location.append((0,i))
    boundary_location.append((width-1,i))

output_matrix = two_dimension_jacobi(metal_box, boundary_location, microVolt, middle_matrix=25)

fig, (ax1,ax2,ax3) = plt.subplots(1,3)
img1 = ax1.imshow(metal_box, cmap="bone")
img2 = ax2.imshow(output_matrix[1], cmap="bone")
img3 = ax3.imshow(output_matrix[0], cmap="bone")
fig.colorbar(img1, ax=[ax1, ax2, ax3], orientation="horizontal")
ax1.title.set_text('Input Matrix')
ax2.title.set_text('Intermediate Matrix')
ax3.title.set_text('Output Matrix')

# Voltage
# jacobi_solution = output_matrix[0]
# coordinates  = np.linspace(-1.0, 1.0, len(jacobi_solution))
# X, Y   = np.meshgrid(coordinates, coordinates)
# dx, dy = np.diff(X, axis=1).item(0), np.diff(Y, axis=0).item(0)
# Ey, Ex = np.gradient(- jacobi_solution)
# Z      = np.sqrt((Ex/dx) ** 2 + (Ey/dy) **2)
# plt.streamplot(X, Y, Ex, Ey, density=1, color='k', linewidth=Z)

plt.show()
