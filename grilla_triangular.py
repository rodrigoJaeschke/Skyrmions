# coding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import funciones_spin as sp
import grilla_rectangular as gr
from scipy.fftpack import fft, fftfreq
from matplotlib.colors import LogNorm


def iniciar_grilla_triangular(nx, ny, j_intercambio, D, cb_periodica=True, d=1, cond_inicial=0, prim_iter=True):
	posiciones = np.zeros((nx * ny , 2))
	for i in range(nx):
		for j in range(ny):
			indice = gr.indice_grilla_rectangular(nx, ny, i, j)
			pos_x = d * i + j * d / 2
			pos_y = j * d * np.sqrt(3) / 2
			posiciones[indice, :] = [pos_x, pos_y]
	listas_vecinos = lista_vecinos_triangular(nx, ny, j_intercambio, D, cb_periodica)
	if prim_iter:
		condicion_inicial = gr.condicion_inicial_aleatoria(nx * ny)
		#condicion_inicial = gr.equilibrio_rectangular(nx, ny, j_intercambio)
	else:
		condicion_inicial = cond_inicial
	
	lista_nodos = []
	for j in range(ny):
		for i in range(nx):
			indice_nodo = gr.indice_grilla_rectangular(nx, ny, i, j)
			nodo = sp.Nodo_spin(condicion_inicial[indice_nodo], 0, listas_vecinos[indice_nodo])
			lista_nodos.append(nodo)
	return lista_nodos, posiciones




def lista_vecinos_triangular(nx, ny, j_intercambio, D, cb_periodica):
	listas_vecinos = []
	for j in range(ny):
		for i in range(nx):
			indices = []
			D_array = []
			for k in range(6):
				D_array = agregar_vecino_triangular(indices, D_array, D, nx, ny, i, j, cb_periodica, k)
			vecinos = np.zeros((len(indices), 5))
			for k in range(len(indices)):
				vecinos[k, 0] = indices[k]
				vecinos[k, 1] = j_intercambio
				Dx = D_array[k][0]
				Dy = D_array[k][1]
				vecinos[k, 2:5] = [Dx, Dy, 0]
			listas_vecinos.append(vecinos)
	return listas_vecinos





def agregar_vecino_triangular(indices, D_array, D_modulo, nx, ny, i, j, cb_periodica, id_vecino):
	D = D_modulo * np.array([1, 0])
	if (id_vecino == 0):
		if not(i==nx - 1):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i + 1, j)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)

		elif cb_periodica:
			indice_vec = gr.indice_grilla_rectangular(nx, ny, 0, j)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)

	if (id_vecino == 1):
		if not(j == ny - 1):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i, j + 1)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)
		elif cb_periodica:
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i, 0)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)

	if (id_vecino == 2):
		if not(j == ny - 1 or i == 0):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i - 1, j + 1)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)
		elif cb_periodica:
			if i == 0 and not(j == ny - 1):
				indice_vec = gr.indice_grilla_rectangular(nx, ny, nx - 1, j + 1)
				indices.append(indice_vec)
			if not(i == 0) and j == ny-1:
				indice_vec = gr.indice_grilla_rectangular(nx, ny, i-1, 0)
				indices.append(indice_vec)
			if i==0 and j==ny - 1:
				indice_vec = gr.indice_grilla_rectangular(nx, ny, i-1, 0)
				indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)

	if (id_vecino == 3):
		if not(i == 0):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i - 1, j)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)
		elif cb_periodica:
			indice_vec = gr.indice_grilla_rectangular(nx, ny, nx - 1, j)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)


	if (id_vecino == 4):
		if not(j == 0):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i, j - 1)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)
		elif cb_periodica:
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i, ny - 1)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)

	if (id_vecino == 5):
		if not(i == nx - 1 or j==0):
			indice_vec = gr.indice_grilla_rectangular(nx, ny, i + 1, j - 1)
			indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)
		elif cb_periodica:
			if i == nx - 1 and not(j ==0):
				indice_vec = gr.indice_grilla_rectangular(nx, ny, 0, j-1)
				indices.append(indice_vec)
			if not(i==nx-1) and j==0:
				indice_vec = gr.indice_grilla_rectangular(nx, ny, i+1, ny - 1)
				indices.append(indice_vec)
			if i==nx-1 and j == 0:
				indice_vec = gr.indice_grilla_rectangular(nx, ny, 0, ny - 1)
				indices.append(indice_vec)
			D = gr.rotar_2d(D, id_vecino * np.pi/3)
			D_array.append(D)


	return D_array




def wending_number_triangular(spines_en_t, nx, ny):
	angulo_solido = 0
	for j in range(ny):
		for i in range(nx):
			area_triangulos = w_n_nodo_triangular(spines_en_t, nx, ny, i, j)
			angulo_solido = angulo_solido + area_triangulos

	return angulo_solido / 4 / np.pi



def w_n_nodo_triangular(spines_en_t, nx, ny, i, j):
	if not(i == nx - 1):
		vec0 = gr.indice_grilla_rectangular(ny, ny, i + 1, j)
	else:
		vec0 = gr.indice_grilla_rectangular(nx, ny, 0, j)


	if not(j == ny - 1):
		vec1 = gr.indice_grilla_rectangular(nx, ny, i, j + 1)
	else:
		vec1 = gr.indice_grilla_rectangular(nx, ny, i, 0)


	if not((i == 0) or (j == ny - 1)):
		vec2 = gr.indice_grilla_rectangular(nx, ny, i - 1, j + 1)
	elif (i == 0) and not(j == ny - 1):
		vec2 = gr.indice_grilla_rectangular(nx, ny, nx - 1, j + 1)
	elif not(i == 0) and (j == ny - 1):
		vec2 = gr.indice_grilla_rectangular(nx, ny, i - 1, 0)
	else:
		vec2 = gr.indice_grilla_rectangular(nx, ny, nx - 1, 0)


	s_nodo = spines_en_t[gr.indice_grilla_rectangular(nx, ny, i, j)]
	s0 = spines_en_t[vec0, :]
	s1 = spines_en_t[vec1, :]
	s2 = spines_en_t[vec2, :]
	area_triangulo_1 = sp.solid_angle_triangle(s_nodo, s0, s1)
	area_triangulo_2 = sp.solid_angle_triangle(s_nodo, s1, s2)

	return area_triangulo_1 + area_triangulo_2




	