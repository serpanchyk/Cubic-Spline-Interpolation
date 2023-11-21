matrix = [
    [1,	9,	0,	0,	0,	0,	0,	0,	0,	0,	0],
    [0,	0,	3,	9,	27,	0,	0,	0,	0,	0,	0],
    [0,	0,	0,	0,	0,	1,	2,	4,	0,	0,	0],
    [0,	0,	0,	0,	0,	0,	0,	0,	4,	16,	64],
    [1,	27,	-1,	0,	0,	0,	0,	0,	0,	0,	0],
    [0,	9,	-1,	0,	0,	0,	0,	0,	0,	0,	0],
    [0,	0,	1,	6,	27,	-1,	0,	0,	0,	0,	0],
    [0,	0,	0,	0,	0,	1,	4,	12,	-1,	0,	0],
    [0,	0,	0,	1,	9,	0,	-1,	0,	0,	0,	0],
    [0,	0,	0,	0,	0,	0,	1,	6,	0,	-1,	0],
    [0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	12]    
]

vector = [1, -4, 1, -3, 0, 0, 0, 0, 0, 0, 0]

def multiplication_by_a_number(number, matrix):
    matrix = [[j*number for j in i] for i in matrix]
    return matrix

def summation_of_matrices(list_of_matrices):
    sum_of_matrices = [[0 for j in range(len(list_of_matrices[0][0]))] for i in range(len(list_of_matrices[0]))]
    #К-ть сумацій матриць, на одне менше від к-ть матриць
    for i in range(len(list_of_matrices)):
        intermediate_sum = []
        #прогонка по рядкам
        for j in range(len(list_of_matrices[i])):
            #сумація по стовпцям
            sum_of_matrices[j] = [x + y for (x, y) in zip(sum_of_matrices[j], list_of_matrices[i][j])]
    return sum_of_matrices

def multiplication_of_matrices(A, B):
    lines_A = len(A)
    lines_B = len(B)
    columns_A = len(A[0])
    columns_B = len(B[0])
    if columns_A == lines_B:
        C = [[0 for line in range(columns_B)] for column in range(lines_A)]
        for i in range(lines_A):
            for j in range(columns_B):
                for n in range(columns_A):
                    C[i][j] += A[i][n] * B[n][j]
        return C
    print('We can`t multiply this matrices')

def permutation_lines(matrix, first_line, second_line):
    i1, i2 = matrix.index(first_line), matrix.index(second_line)
    matrix[i1], matrix[i2] = matrix[i2].copy(), matrix[i1].copy()
    return matrix

def transpose_matrix(matrix): 
    if type(matrix[0]) != list:
        transpose_matrix = [[i] for i in matrix] 
        return transpose_matrix
    transpose_matrix = [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]
    return transpose_matrix

def get_determinant(matrix):
    lenght = len(matrix)
    i1, i2 = 1, 1
    determinant = 1
    permutation = 0
    #Перевірка на нульове значення першого елемента та усунення цього
    while matrix[0][0] == 0:
        if i1 == lenght:
            if i2 == lenght:
                break
            else:
                matrix[0] = [x+y for x, y in zip(matrix[0], matrix[i2])]
                i2 += 1  
        else: 
            matrix = permutation_lines(matrix.copy(), matrix[0], matrix[i1])
            permutation += 1
            i1 += 1 

    for i in range(1, lenght):
        for j in range(i, lenght):
            if matrix[i-1][i-1] == 0:
                for line in range(i, lenght):
                    if matrix[line][i-1] != 0:
                        matrix = permutation_lines(matrix.copy(), matrix[i-1], matrix[line]) 
                        permutation += 1
                        matrix[j] = [x-y*(matrix[j][i-1]/matrix[i-1][i-1]) for x, y in zip(matrix[j], matrix[i-1])] 
            else:  
                matrix[j] = [x-y*(matrix[j][i-1]/matrix[i-1][i-1]) for x, y in zip(matrix[j], matrix[i-1])]
    for i in range(lenght):
        determinant *= matrix[i][i]
    if permutation != 0:
        for i in range(0, permutation):
            determinant *= -1
    return determinant


def get_minor(matrix, y, x):
    minor = [line[:x] + line[x+1:] for line in (matrix[:y]+matrix[y+1:])]
    return minor


def algebraic_complement(matrix, y, x):  
    algebraic_complement = (-1) ** (y+x) * get_determinant(get_minor(matrix.copy(), y, x))
    return algebraic_complement
    

def get_cofactor_matrix(matrix):
    lenght = len(matrix)
    cofactor_matrix = [i*0 for i in range(lenght)]
    for i in range(lenght):
        cofactor_matrix[i] = [j*0 for j in range(lenght)]
    for y in range(lenght):
        for x in range(lenght):
            cofactor_matrix[y][x] = algebraic_complement(matrix.copy(), y, x)
    return cofactor_matrix


def get_inverse_matrix(matrix):
    det = get_determinant(matrix.copy())
    if det == 0:
        print('We haven`t found inverse matrix')
    else:
        inverse_matrix = multiplication_by_a_number(get_determinant(matrix.copy()) ** -1, transpose_matrix(get_cofactor_matrix(matrix.copy())))
    return inverse_matrix


def bipolarization(vector):
    for i in range(len(vector)):
        if vector[i] == 0:
            vector[i] = -1
    return vector


def unipolarization(vector):
    for i in range(len(vector)):
        if vector[i] == -1:
            vector[i] = 0
    return vector

        
def gauss(matrix, vector):
    epsilon = 0.0000000001
    lenght = len(matrix)
    coefficients = [0 for i in range(lenght)]
    transpose_matrix = [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix))]
    for i in range(lenght):
        if len(matrix[i]) != lenght:
            print('Matrycja ne kvadratna!')
            return coefficients   
    for i in range(lenght):
        if sum(matrix[i]) == 0:
            print('Matrycja ne kvadratna!')
            return coefficients
        elif sum(transpose_matrix[i]) == 0:
            print('Matrycja ne kvadratna!')
            return coefficients   
    for i in range(lenght):
        if abs(matrix[0][0]) > epsilon:
            break
        else:
            matrix[0], matrix[i] = matrix[i][:], matrix[0][:]
    for i in range(lenght):
        if matrix[i][i] == 0:
            for q in range(i, lenght):
                if abs(matrix[q][i]) > epsilon:
                    matrix[q], matrix[i] = matrix[i][:], matrix[q][:]
                    vector[q], vector[i] = vector[i], vector[q]
                    break
        number = matrix[i][i] ** -1        
        matrix[i] = [x*number for x in matrix[i]]  
        vector[i] *= number 
        for k in range(i+1, lenght):
            number = -matrix[k][i]
            matrix[k] = [x+y*number for x, y in zip(matrix[k], matrix[i])]
            vector[k] += vector[i]*number 
    for j in range(lenght-1, -1, -1):
        for i in range(j-1, -1, -1):
            number = -matrix[i][j]
            matrix[i] = [x+y*number for x, y in zip(matrix[i], matrix[j])]
            vector[i] += vector[j]*number
    coefficients = vector.copy()
    return coefficients