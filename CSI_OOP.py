from operator import __add__, __sub__, __mul__
from math_matrix import get_inverse_matrix, multiplication_of_matrices
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from math import sin, cos, pi

import numpy as np
from scipy.interpolate import CubicSpline

x = 0
y = 1
z = 2

# клас кубічних сплайнів
class CSI:
    def __init__(self, Pk, P1_d, Pn_d, dimensions):
        self.number_of_points = len(Pk)
        self.number_of_segments = len(Pk) - 1
        self.Pk = Pk
        self.Pk1_d = P1_d
        self.Pkn_d = Pn_d
        self.dimensions = dimensions

    def interpolation(self, boundary_conditions='fixed', colour='b', tangent_vectors=False, tangent_vectors_colour='k', dots=False, dots_colour='k'):
        # Задавання параметра хордовою інтерполяцією
        t = []
        if self.dimensions == 2:
            for i in range(self.number_of_segments):
                t.append(((self.Pk[i+1][x] - self.Pk[i][x]) ** 2 + (self.Pk[i+1][y] - self.Pk[i][y]) ** 2) ** 0.5)
        elif self.dimensions == 3:
            for i in range(self.number_of_segments):
                t.append(((self.Pk[i+1][x] - self.Pk[i][x]) ** 2 + (self.Pk[i+1][y] - self.Pk[i][y]) ** 2 + (self.Pk[i+1][z] - self.Pk[i][z]) ** 2) ** 0.5)    
        
        # Задання матриць для пошуку дотичних векторів до сегментів
        M = [[0 for j in range(self.number_of_points)] for i in range(self.number_of_points-2)]
        for i in range(self.number_of_points-2):
            M[i][i] = t[i+1]
            M[i][i+1] = 2 * (t[i+1] + t[i])
            M[i][i+2] = t[i]
        a, b = [0 for i in range(self.number_of_points)], [0 for i in range(self.number_of_points)]     
        if boundary_conditions == 'fixed':
            a[0], b[len(b)-1] = 1, 1
            M.insert(len(b)-1, b)
        elif boundary_conditions == 'weak':
            a[0], a[1], b[len(b)-2], b[len(b)-1] = 1, 1/2, 2, 4
            M.insert(len(b)-1, b)
        elif boundary_conditions == 'cyclic':
            a[0], a[1], a[len(a)-2] = 2*(1 + t[self.number_of_segments-1]/t[0]), t[self.number_of_segments-1]/t[0], 1
        elif boundary_conditions == 'acyclic':
            a[0], a[1], a[len(a)-2] = 2*(1 + t[self.number_of_segments-1]/t[0]), t[self.number_of_segments-1]/t[0], -1    
        M.insert(0, a)
        

        R = []
        if boundary_conditions == 'fixed':
            R.append(self.Pk1_d)
        elif boundary_conditions == 'weak':
            R.append([3/(2*t[0])*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))])
        elif boundary_conditions == 'cyclic':
            first_subtraction = [3*t[self.number_of_segments-1]/t[0]**2*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))]
            second_subtraction = [3/t[self.number_of_segments-1]*j for j in list(map(__sub__, self.Pk[self.number_of_segments-1], self.Pk[self.number_of_segments]))]
            R.append(list(map(__sub__, first_subtraction, second_subtraction)))
        elif boundary_conditions == 'acyclic':
            first_subtraction = [3*t[self.number_of_segments-1]/t[0]**2*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))]
            second_subtraction = [3/t[self.number_of_segments-1]*j for j in list(map(__sub__, self.Pk[self.number_of_segments-1], self.Pk[self.number_of_segments]))]
            R.append(list(map(__add__, first_subtraction, second_subtraction)))               
        
        for i in range(self.number_of_points-2):
            first_subtraction = [t[i]**2 * j for j in list(map(__sub__, self.Pk[i+2], self.Pk[i+1]))]
            second_subtraction = [t[i+1] ** 2 * j for j in list(map(__sub__, self.Pk[i+1], self.Pk[i]))]
            summation = list(map(__add__, first_subtraction, second_subtraction))
            division = 3 / (t[i]*t[i+1])
            multiplication = [division * j for j in summation]
            R.append(multiplication)
        
        if boundary_conditions == 'fixed':
            R.append(self.Pkn_d)
        elif boundary_conditions == 'weak':
            R.append([6/t[self.number_of_segments-1]*j for j in list(map(__sub__, self.Pk[self.number_of_segments], self.Pk[self.number_of_segments-1]))])
        
        M_inverse = get_inverse_matrix(M)

        # Добуток матриць - дотичні вектори до сегментів
        P_d = multiplication_of_matrices(M_inverse, R)
        if boundary_conditions == 'cyclic':
            P_d.append(P_d[0][:])
        elif boundary_conditions == 'acyclic':
            P_d.append([i*(-1) for i in P_d[0][:]])

        # Знаходження точок на сегменті
        step = 0.01
        P = []
        for i in range(self.number_of_segments):
            r = 0
            while r <= 1:
                F = [[2*r**3-3*r**2+1, -2*r**3+3*r**2, r*(r**2-2*r+1)*t[i], r*(r**2-r)*t[i]]]
                G = [self.Pk[i], self.Pk[i+1], P_d[i], P_d[i+1]]
                P.append(multiplication_of_matrices(F, G)[0])
                r += step   

        if self.dimensions == 2:
            if tangent_vectors == True:
                for i in range(len(P_d)):
                    ax.arrow(self.Pk[i][x], self.Pk[i][y], P_d[i][x], P_d[i][y], head_width=0.05, head_length=0.1, fc=tangent_vectors_colour, ec=tangent_vectors_colour)
            Px = []
            Py = []
            for i in P:
                Px.append(i[x])
                Py.append(i[y])
            ax.plot(Px, Py, f'{colour}')
            if dots == True:
                for i in range(len(self.Pk)):
                    ax.plot(self.Pk[i][x], self.Pk[i][y], dots_colour)

        elif self.dimensions == 3:
            if tangent_vectors == True:
                for i in range(len(P_d)):
                    ax.quiver(self.Pk[i][x], self.Pk[i][y], self.Pk[i][z], P_d[i][x], P_d[i][y], P_d[i][z], color=tangent_vectors_colour)
            Px = []
            Py = []
            Pz = []
            for i in P:
                Px.append(i[x])
                Py.append(i[y])
                Pz.append(i[z])
            ax.plot3D(Px, Py, Pz, f'{colour}')
            if dots == True:
                for i in range(len(self.Pk)):
                    ax.plot(self.Pk[i][x], self.Pk[i][y], self.Pk[i][z], dots_colour)
        
        return(P)


    def normalise_interpolation(self, boundary_conditions='fixed', colour='b', tangent_vectors=False, tangent_vectors_colour='k', dots=False, dots_colour='g'):
        # Задання матриць для пошуку дотичних векторів до сегментів
        M = [[0 for j in range(self.number_of_points)] for i in range(self.number_of_points-2)]
        for i in range(self.number_of_points-2):
            M[i][i] = 1
            M[i][i+1] = 4
            M[i][i+2] = 1
        a, b = [0 for i in range(self.number_of_points)], [0 for i in range(self.number_of_points)]     
        if boundary_conditions == 'fixed':
            a[0], b[len(b)-1] = 1, 1
            M.insert(len(b)-1, b)
        elif boundary_conditions == 'weak':
            a[0], a[1], b[len(b)-2], b[len(b)-1] = 1, 1/2, 2, 4
            M.insert(len(b)-1, b)
        elif boundary_conditions == 'cyclic':
            a[0], a[1], a[len(a)-2] = 2*(1 + 1), 1, 1
        elif boundary_conditions == 'acyclic':
            a[0], a[1], a[len(a)-2] = 2*(1 + 1), 1, -1    
        M.insert(0, a)

        R = []
        if boundary_conditions == 'fixed':
            R.append(self.Pk1_d)
        elif boundary_conditions == 'weak':
            R.append([3/(2*1)*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))])
        elif boundary_conditions == 'cyclic':
            first_subtraction = [3*1**2*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))]
            second_subtraction = [3/1*j for j in list(map(__sub__, self.Pk[self.number_of_segments-1], self.Pk[self.number_of_segments]))]
            R.append(list(map(__sub__, first_subtraction, second_subtraction)))
        elif boundary_conditions == 'acyclic':
            first_subtraction = [3*1**2*j for j in list(map(__sub__, self.Pk[1], self.Pk[0]))]
            second_subtraction = [3/1*j for j in list(map(__sub__, self.Pk[self.number_of_segments-1], self.Pk[self.number_of_segments]))]
            R.append(list(map(__add__, first_subtraction, second_subtraction)))               
        
        for i in range(self.number_of_points-2):
            first_subtraction = list(map(__sub__, self.Pk[i+2], self.Pk[i+1]))
            second_subtraction = list(map(__sub__, self.Pk[i+1], self.Pk[i]))
            summation = list(map(__add__, first_subtraction, second_subtraction))
            multiplication = [3 * j for j in summation]
            R.append(multiplication)
        
        if boundary_conditions == 'fixed':
            R.append(self.Pkn_d)
        elif boundary_conditions == 'weak':
            R.append([6/1*j for j in list(map(__sub__, self.Pk[self.number_of_segments], self.Pk[self.number_of_segments-1]))]) 
        
        M_inverse = get_inverse_matrix(M)
        
        # Добуток матриць - дотичні вектори до сегментів
        P_d = multiplication_of_matrices(M_inverse, R)
        if boundary_conditions == 'cyclic':
            P_d.append(P_d[0][:])
        elif boundary_conditions == 'acyclic':
            P_d.append([i*(-1) for i in P_d[0][:]])

        # Знаходження точок на сегменті
        step = 0.01
        P = []
        for i in range(self.number_of_segments):
            t = 0
            while t <= 1:
                T = [[t ** j for j in range(3, -1, -1)]]
                N = [
                    [2, -2, 1, 1],
                    [-3, 3, -2, -1],
                    [0, 0, 1, 0],
                    [1, 0, 0, 0]
                    ]
                F = multiplication_of_matrices(T, N)
                G = [self.Pk[i], self.Pk[i+1], P_d[i], P_d[i+1]]
                P.append(multiplication_of_matrices(F, G)[0])
                t += step
        
        if self.dimensions == 2:
            if tangent_vectors == True:
                for i in range(len(P_d)):
                    ax.arrow(self.Pk[i][x], self.Pk[i][y], P_d[i][x], P_d[i][y], head_width=0.05, head_length=0.1, fc=tangent_vectors_colour, ec=tangent_vectors_colour)
            Px = []
            Py = []
            for i in P:
                Px.append(i[x])
                Py.append(i[y])
            ax.plot(Px, Py, f'{colour}')
            if dots == True:
                for i in range(len(self.Pk)):
                    ax.plot(self.Pk[i][x], self.Pk[i][y], dots_colour)

        elif self.dimensions == 3:
            if tangent_vectors == True:
                for i in range(len(P_d)):
                    ax.quiver(self.Pk[i][x], self.Pk[i][y], self.Pk[i][z], P_d[i][x], P_d[i][y], P_d[i][z], color=tangent_vectors_colour)
            Px = []
            Py = []
            Pz = []
            for i in P:
                Px.append(i[x])
                Py.append(i[y])
                Pz.append(i[z])
            ax.plot3D(Px, Py, Pz, f'{colour}')
            if dots == True:
                for i in range(len(self.Pk)):
                    ax.plot(self.Pk[i][x], self.Pk[i][y], self.Pk[i][z], dots_colour)
        
        return(P)

# Ця зміна відповідає за вибір між двовимірним і тривимірним зображенням 
dimensions = 2
if dimensions == 2:
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.axis('equal')
elif dimensions == 3:
    ax = plt.axes(projection='3d')
    ax.set_zlabel("z", fontsize=14)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-8, 8)
    ax.set_zlim(-8, 8)
ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("y", fontsize=14)
# ax.grid(which="major", linewidth=1.2)
# ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.tick_params(which='major', length=10, width=2)
# ax.tick_params(which='minor', length=5, width=1)


P1k = [[0, 4.1], [1, 2.4], [2, 3], [3, 4.3], [4, 3.6], [5, 5.2], [6, 5.9]]
P1_d = [0, 0]
Pn_d = [0, 0]
spline = CSI(P1k, P1_d, Pn_d, dimensions)
spline.interpolation(colour='b', dots=True, dots_colour='bo', boundary_conditions='fixed')

x = np.arange(0, 6)



# P1k = [[0, 4.1], [1, 2.4], [2, 3], [3, 4.3], [4, 3.6], [5, 5.2], [6, 5.9]]
# P1_d = [0, -1]
# Pn_d = [1, 0]
# spline = CSI(P1k, P1_d, Pn_d, dimensions)
# spline.normalise_interpolation(colour='b', dots=True, dots_colour='bo', boundary_conditions='fixed', tangent_vectors=True, tangent_vectors_colour='k')
# spline.normalise_interpolation(colour='g', dots=True, dots_colour='go', boundary_conditions='weak', tangent_vectors=True, tangent_vectors_colour='k')
# spline.normalise_interpolation(colour='y', dots=True, dots_colour='yo', boundary_conditions='cyclic', tangent_vectors=True, tangent_vectors_colour='k')
# spline.interpolation(colour='b', dots=True, dots_colour='bo', boundary_conditions='fixed', tangent_vectors=True, tangent_vectors_colour='grey')

# Модель келиха
# P1k = [[0, 0, 0], [3, 0, 0], [3, 0, 0.3], [0.7, 0, 1.5], [0.7, 0, 7.5], [3.7, 0, 10], [4, 0, 15], [3.3, 0, 10], [0, 0, 8.2]]
# P1_d = [0, 0, 0]
# Pn_d = [0, 0, 0]
# a = CSI(P1k, P1_d, Pn_d, dimensions)
# P = a.normalise_interpolation(colour='b')
# ax.set_xlim(-8, 8)
# ax.set_ylim(-8, 8)
# ax.set_zlim(0, 16)

# pieces = 10
# step = 2*pi / pieces
# xi = 0
# for i in range(pieces):  
#     Rz = [
#         [cos(xi), -sin(xi), 0],
#         [sin(xi), cos(xi), 0],
#         [0, 0, 1]
#         ]
#     xi += step    
#     Px = []
#     Py = []
#     Pz = []
#     for i in P:
#         P_r = multiplication_of_matrices([i], Rz)[0]
#         Px.append(P_r[x])
#         Py.append(P_r[y])
#         Pz.append(P_r[z])
#     ax.plot3D(Px, Py, Pz, 'b')

# step = 50
# for i in range(0, len(Px), step):
#     r = P[i][x]
#     t = 0
#     x_circle = []
#     y_circle = []
#     while t < 2*pi:
#         x_circle.append(r*cos(t))
#         y_circle.append(r*sin(t))
#         t += 0.01
#     z_circle = [P[i][z] for j in range(len(x_circle))]
#     ax.plot3D(x_circle, y_circle, z_circle, 'b')




# Інтерполяція гвинтової лінії
# ax.set_xlim(-4, 4)
# ax.set_ylim(-4, 4)
# ax.set_zlim(0, 16)
# P1_d = [0, 3.5, 0]
# Pn_d = [0, 3.5, 0]
# Px = []
# Py = []
# Pz = []
# a = 2
# b = 1/2
# step1 = 0.01
# h = 5*pi/b
# for t in range(int(h/step1)+1):
#     t /= (1/step1)
#     Px.append(a*cos(t))
#     Py.append(a*sin(t))
#     Pz.append(b*t)
# ax.plot3D(Px, Py, Pz, 'r')

# P2k = []
# step2 = pi/4
# for i in range(0, int(h/step1)+1, int(step2/step1/b)):
#     P2k.append([Px[i], Py[i], Pz[i]])
# spline = CSI(P2k, P1_d, Pn_d, dimensions)
# spline.normalise_interpolation(colour='b', dots=True, dots_colour='bo', boundary_conditions='fixed')




# Інтерполяція тороїдальної спіралі
# Px = []
# Py = []
# Pz = []
# R = 6
# r = 1
# n = 5
# t = 0
# step = 0.01
# while t < 2*pi:
#     t += step
#     Px.append((R+r*cos(n*t))*cos(t))
#     Py.append((R+r*cos(n*t))*sin(t))
#     Pz.append(pi*r*sin(n*t))
# ax.plot3D(Px, Py, Pz, 'r')

# P1k = []
# dots = 19
# k = 0
# for i in range(0, len(Px)-int(len(Px)/dots), int(len(Px)/dots)):
#     k += 1
#     P1k.append([Px[i], Py[i], Pz[i]])
# P1k.append([Px[0], Py[0], Pz[0]])
# P1_d = [0, 0, 0]
# Pn_d = [0, 0, 0]

# spline = CSI(P1k, P1_d, Pn_d, dimensions)
# spline.interpolation(colour='g', dots=True, dots_colour='go', boundary_conditions='fixed')

plt.show()