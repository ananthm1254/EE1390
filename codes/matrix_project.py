# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import subprocess

def tangent(p):
	n0 = np.matmul(p,V)
	k0 = 1
	omat = np.array([[0,-1],[1,0]])
	m0 = np.matmul(omat,n0.T)
	len = 100
	tang_x = np.zeros((2,len))
	lam = np.linspace(-10,10,len)
	for i in range(0,100):
		temp = p + lam[i]*(m0.T)
		tang_x[:,i] = temp.T
	return tang_x

def line_intersect(n1,n2):
    N = np.vstack((n1, n2))
    p = np.zeros(2)
    p[0] = 1
    p[1] = 1
    return np.matmul(np.linalg.inv(N),p)

#features of the eclipse V and 
V = np.array([[1.0/9, 0],[0, 1.0/(5)]])
F = -1

p1 = np.array([2.0,5.0/3])
p2 = np.array([2.0,-5.0/3])
p3 = np.array([-2.0,5.0/3])
p4 = np.array([-2.0,-5.0/3])

plt.plot(p1[0],p1[1],'x')
plt.plot(p2[0],p2[1],'x')
plt.plot(p3[0],p3[1],'x')
plt.plot(p4[0],p4[1],'x')

tang_1 = tangent(p1)
tang_2 = tangent(p2)
tang_3 = tangent(p3)
tang_4 = tangent(p4)    

plt.plot(tang_1[0,:],tang_1[1,:], label = '$6x + 9y = 27$')
plt.plot(tang_2[0,:],tang_2[1,:], label = '$6x - 9y = 27$')
plt.plot(tang_3[0,:],tang_3[1,:], label = '$-6x - 9y = 27$')
plt.plot(tang_4[0,:],tang_4[1,:], label = '$-6x + 9y = 27$')

n_3 = np.matmul(p3,V)
n_4 = np.matmul(p4,V)
n_1 = np.matmul(p1,V)
n_2 = np.matmul(p2,V)

A = line_intersect(n_1,n_2)
B = line_intersect(n_2,n_4)
C = line_intersect(n_3,n_4)
D = line_intersect(n_3,n_1)

print('\nPoints of intersection of tangents:')
print(A)
print(B)
print(C)
print(D)

plt.plot(A[0],A[1],'o')
plt.text(A[0]*(1+0.1),A[1],'$A$')
plt.plot(B[0],B[1],'o')
plt.text(B[0],B[1]*(1+0.1),'$B$')
plt.plot(C[0],C[1],'o')
plt.text(C[0]*(1+0.1),C[1],'C')
plt.plot(D[0],D[1],'o')
plt.text(D[0],D[1]*(1+0.1),'D')

area = (2*A[0]*D[1])
print('\n')
print('Area of the quadrilateral required =')
print(area)

len = 100000
t = np.linspace(0,2*np.pi, 100000)
eclipse = np.zeros((2,len))
A = np.array([[3,0],[0,np.sqrt(5)]])
for i in range(0,len):
    temp = np.matmul(np.array([np.cos(t[i]), np.sin(t[i])]),A)
    eclipse[:,i] = temp.T 
plt.plot(eclipse[0,:],eclipse[1,:], label='$x^2/9 + y^2/5 = 1$')
plt.grid()
plt.legend(loc='upper right')
plt.show()
    
