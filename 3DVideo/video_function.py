import numpy as np

# Variables for testing delta V 
deciVolt  = 0.1
centiVolt = 0.01
milliVolt = 0.001
microVolt = 0.000001
nanoVolt  = 0.000000001
picoVolt  = 0.000000000001
femtoVolt = 0.000000000000001
attoVolt  = 0.000000000000000001
zeptoVolt = 0.000000000000000000001
yoctoVolt = 0.000000000000000000000001

def two_dimension_jacobi(schematic, boundary_location, error_treshold):
    height, width = schematic.shape

    # As input_matrix is updated, the result is stored in the output_matrix.
    input_matrix  = schematic.copy()
    output_matrix = schematic.copy()

    # The number of iterations to reach acceptable delta_v is tracked.
    iterations = 0
    delta_v = 0

    all_matrices = []

    while True:
        # Loop left to right from top to bottom
        all_matrices.append([input_matrix])
        for i in range(1, width -1):
            for j in range(1, height -1):
                if (i,j) in boundary_location:
                    pass
                else:  # New voltage value = average of left neighbour, right neighbour, bottom neighbour, top neighbour values
                    output_matrix[i,j] = (input_matrix[i-1][j] + input_matrix[i+1][j] + input_matrix[i][j-1] + input_matrix[i][j+1])/4
                    delta_v += abs(np.subtract(input_matrix[i,j], output_matrix[i,j]))

        iterations += 1

        if delta_v < error_treshold:
            break
        else:
            delta_v = 0  # Restart counting delta_v for the next iteration 
            input_matrix = output_matrix.copy() # Deep copy input_matrix to output_matrix

    return output_matrix, iterations, delta_v, all_matrices


def two_dimension_gauss_seidel(schematic, boundary_location, error_treshold):
    # Get the dimensions of the input matrix to loop through. 
    height, width = schematic.shape

    # As input_matrix is updated, the result is stored in the input_matrix (not the output_matrix)
    # A copy of the input matrix is kept to make calculating delta_v easy and intuitive
    input_matrix        = schematic.copy()
    input_matrix_copy   = input_matrix.copy()

    # The number of iterations to reach acceptable delta_v is tracked.
    iterations = 0
    delta_v = 0

    while True:
        # Loop left to right from top to bottom
        for i in range(1, width -1):
            for j in range(1, height -1):
                if (i,j) in boundary_location:
                    pass
                else:  # New voltage value = average of left neighbour, right neighbour, bottom neighbour, top neighbour values
                    input_matrix[i,j] = (input_matrix[i-1][j] + input_matrix[i+1][j] + input_matrix[i][j-1] + input_matrix[i][j+1])/4
                    delta_v += abs(np.subtract(input_matrix[i,j], input_matrix_copy[i,j]))

        iterations += 1
        if delta_v < error_treshold:
            break
        else:
            delta_v = 0  # Restart counting delta_v for the next iteration 
            input_matrix_copy = input_matrix.copy() # Deep copy input_matrix to output_matrix

    return input_matrix, iterations, delta_v


def two_dimension_sor(schematic, boundary_location, alpha, error_treshold):
    # Get the dimensions of the input matrix to loop through. 
    height, width = schematic.shape
    # alpha = 2/(1+(np.pi/width))

    # As the input_matrix is updated, the result is stored in the input_matrix (Gauss Seidel solution).
    # A new array is made for the output_matrix as it does not equal the Gauss Seidel matrix
    input_matrix        = schematic.copy()
    gauss_seidel        = schematic.copy()
    output_matrix       = schematic.copy()

    # The number of iterations to reach acceptable delta_v is tracked.
    iterations = 0
    delta_v = 0

    while True:
        # Loop left to right from top to bottom
        for i in range(1, width -1):
            for j in range(1, height -1):
                if (i,j) in boundary_location:
                    pass
                else:  # New voltage value = average of left neighbour, right neighbour, bottom neighbour, top neighbour values
                    gauss_seidel[i,j] = (input_matrix[i-1][j] + input_matrix[i+1][j] + input_matrix[i][j-1] + input_matrix[i][j+1])/4
                    delta_v += np.subtract(gauss_seidel[i,j], input_matrix[i,j])
                    output_matrix[i,j] = (alpha * delta_v) + input_matrix[i,j]

        iterations += 1
        if delta_v < error_treshold:
            break
        else:
            delta_v = 0  # Restart counting delta_v for the next iteration 
            input_matrix = gauss_seidel.copy()

    return output_matrix, iterations, delta_v

