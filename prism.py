from matplotlib import pyplot as plt
import numpy as np
from numerical_methods import *

width, height  = 30, 30
prism = np.zeros((width, height), dtype=float)
boundary_location = []
for i in range(int(width * 0.4), int(width * 0.6)):
	for j in range(int(width * 0.4), int(width * 0.6)): 
		prism[i][j] = 1
		boundary_location.append((i,j))



output_matrix = two_dimension_jacobi(prism, boundary_location, microVolt, middle_matrix=25)

fig, (ax1,ax2,ax3) = plt.subplots(1,3)
first = ax1.imshow(prism, cmap="bone")
first = ax2.imshow(output_matrix[1], cmap="bone")
first = ax3.imshow(output_matrix[0], cmap="bone")
fig.colorbar(first, ax=[ax1, ax2, ax3], orientation="horizontal")


## Show Voltage
#jacobi_solution = output_matrix[0]
#coordinates  = np.linspace(0, 1.0, len(jacobi_solution))
#X, Y   = np.meshgrid(coordinates, coordinates)
#dx, dy = np.diff(X, axis=1).item(0), np.diff(Y, axis=0).item(0)
#Ey, Ex = np.gradient(- jacobi_solution)
#Z      = np.sqrt((Ex/dx) ** 2 + (Ey/dy) **2)
##np.savetxt("Z_prism.csv", Z, delimiter=",")
#plt.streamplot(X, Y, Ex, Ey, density=1, color='k', linewidth=Z/2)

plt.show()
