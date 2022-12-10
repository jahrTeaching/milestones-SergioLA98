# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 23:18:00 2022

@author: serg_
"""

from numpy import zeros, size, array, dot
from numpy.linalg import inv, norm

def Jac(F, U):
	N = size(U)
	J= zeros([N,N])
	t = 1e-10

	for i in range(N):
		xj = zeros(N)
		xj[i] = t
		J[:,i] = (F(U + xj) - F(U - xj))/(2*t)
	return J  

def factLU(A):

	N = size(A,1)
	U = zeros([N,N])
	L = zeros([N,N])

	U[0,:] = A[0,:]
	for k in range(0,N):
		L[k,k] = 1

	L[1:N,0] = A[1:N,0]/U[0,0]


	for k in range(1,N):

		for j in range(k,N):
			U[k,j] = A[k,j] - dot(L[k,0:k], U[0:k,j])

		for i in range(k+1,N):
			L[i,k] =(A[i,k] - dot(U[0:k,k], L[i,0:k])) / (U[k,k])

	return [L@U, L, U]


def solveLU(M,b):

	N=size(b)
	y=zeros(N)
	x=zeros(N)

	[A,L,U] = factLU(M)
	y[0] = b[0]

	for i in range(0,N):
		y[i] = b[i] - dot(A[i,0:i], y[0:i])
		

	x[N-1] = y[N-1]/A[N-1,N-1]

	for i in range(N-2,-1,-1):
		x[i] = (y[i] - dot(A[i, i+1:N+1], x[i+1:N+1])) / A[i,i]
		
	return x


def Inv(A):

	N = size(A,1)

	B = zeros([N,N])

	for i in range(0,N):
		one = zeros(N)
		one[i] = 1

		B[:,i] = solveLU(A, one)

	return B


def newton(func, U_0):
	N = size(U_0) 
	U = zeros(N)
	U1 = U_0
	error = 1
	stop = 1e-8
	iteration = 0

	while error > stop and iteration < 1000:
		U = U1 - dot(Inv(Jac(func, U1)),func(U1))
		error = norm(U - U1)
		U1 = U
		iteration = iteration +1
	return U