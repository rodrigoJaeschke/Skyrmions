# coding: utf-8
from __future__ import division
import numpy as np
import funciones_spin as sp
from scipy.fftpack import fft, fftfreq, fft2, fftshift
import random as rand

def iniciar_grilla_rectangular(nx, ny, j, modo='aleatorio'):
	'''
	iniciar_grilla_rectangular(int nx, int ny, float j):
	La funcion recibe  dos dimensiones nx y ny, y genera
	la lista de nodos inicial de una grilla rectangular
	de elementos Nodo_spin con una cte de intercambio "j".

	Ademas de la lista de nodos, retorna una matriz
	con las posiciones en formato [pos_x, pos_y] de
	dimensiones (n_nodos, 2). 
	Se asume una separacion de 4 entre los nodos,
	y ausencia de campo externo.

	return lista_nodos, posiciones
	'''
	H_ext = 0
	separacion = 1
	x = np.linspace(0, (nx - 1) * separacion, nx)
	y = np.linspace(0, (ny - 1) * separacion, ny)
	X, Y = np.meshgrid(x, y)
	posiciones = np.zeros((nx * ny, 2))
	listas_vecinos = lista_vecinos_rectangular(nx, ny, j)

	if (modo == 'peq_osc'):
		condicion_inicial= equilibrio_rectangular(nx, ny, j)
	else:
		condicion_inicial = condicion_inicial_aleatoria(nx * ny)
	lista_nodos = []
	for k in range(nx):
		for l in range(ny):
			i = indice_grilla_rectangular(nx, k, l)
			posiciones[i, :] = [X[k, l], Y[k, l]]
			nodo_nuevo = sp.Nodo_spin(condicion_inicial[i], H_ext, listas_vecinos[i])
			lista_nodos.append(nodo_nuevo)
	return lista_nodos, posiciones


def perturbacion(epsilon):
	ds_x = (rand.random() - 0.5) * epsilon 
	ds_y = (rand.random() - 0.5) * epsilon
	ds_z = (rand.random() - 0.5) * epsilon
	ds = np.array([ds_x, ds_y, ds_z])
	return ds


def equilibrio_rectangular(nx, ny, j, s_eq=np.array([3, 1, 3]), epsilon=0.001):
	cond_inicial = np.zeros((nx * ny, 3))
	s_eq = sp.normalizar(s_eq)
	for k in range(nx):
		for l in range(ny):
			i = indice_grilla_rectangular(nx, k, l)
			ds = perturbacion(epsilon)
			if (j > 0):
				s = sp.normalizar(s_eq * (-1) ** (k + l) + ds)
			else:
				s = sp.normalizar(s_eq + ds)
			cond_inicial[i] = s
	return cond_inicial




def lista_vecinos_rectangular(nx, ny, j_intercambio):
	'''
	lista_vecinos_rectangular(int nx, int ny, float j_intercambio)
	La funcion recibe el taman~o nx y ny de una grilla rectangular
	y la constante j_intercambio entre los spines.
	Retorna una lista 'listas_vecinos' que contiene como elementos,
	a las listas de vecinos de cada Nodo_spin de la grilla rectangular.
	Si un Nodo esta en el borde, se asume condicion periodica y se
	enlaza con el nodo del otro extremo, en la dimension correspondiente.


	'''
	listas_vecinos = []
	for i in range(nx):
	    for j in range(ny):
	        indices = []
	        if not(j == 0):
	            indices.append(indice_grilla_rectangular(nx, i, j - 1))
	        else:
	            indices.append(indice_grilla_rectangular(nx, i, -1))
	            
	        if not(i == nx - 1):
	            indices.append(indice_grilla_rectangular(nx, i + 1, j))
	        else:
	            indices.append(indice_grilla_rectangular(nx, 0, j))
	            
	        if not(j == ny - 1):
	            indices.append(indice_grilla_rectangular(nx, i, j + 1))
	        else:
	            indices.append(indice_grilla_rectangular(nx, i, 0))
	            
	        if not(i == 0):
	            indices.append(indice_grilla_rectangular(nx, i - 1, j))
	        else:
	            indices.append(indice_grilla_rectangular(nx, -1, j))
	            
	        vecinos = np.zeros((len(indices), 2))
	        for k in range(len(indices)):
	            vecinos[k, :] = [indices[k], j_intercambio]
	        listas_vecinos.append(vecinos)
	return listas_vecinos


def indice_grilla_rectangular(nx, i, j):
	'''
	Recibe un taman~o nx (horizontal) de una grilla
	rectangular de nodos, asi como los indices
	(i, j) que definen las coordenadas de un determinado
	Nodo_spin en la grlla.

	La función retorna el índice del correspondiente nodo
	(i, j), con el que se enumerarian los nodos en un arreglo
	unidimensional.
	'''
	return nx * j + i


def condicion_inicial_aleatoria(n):
	'''
	Recibe un entero n, correspondiente al número de nodos
	que se quieren inicialiar. 
	La funcion retorna un arreglo bidimensional
	de dimension (n, 3), que contiene las condiciones
	iniciales x, y, z de 'n' nodos inicializadas de manera
	aleatoria.
	'''
	cond_inicial = np.zeros((n, 3))
	for i in range(n):
	    theta = rand.random() * np.pi
	    phi = rand.random() * 2 * np.pi
	    spin = sp.esf_to_cart(1, theta, phi)
	    cond_inicial[i] = spin
	return cond_inicial


def generar_grilla_peq_osc(lista_nodos0, tolerancia_eq, amp_med_p_o):
	'''
	'''
	t_array = np.linspace(0, 40, 400)
	spin_actual = sp.extraer_spines([lista_nodos0])[:, 0, :]
	dif = tolerancia_eq + 1
	lista_nodos_actual = lista_nodos0
	cont = 0
	grilla_peq_osc = []
	while (dif > tolerancia_eq):
	    grilla_nodos = sp.simular_spines(t_array, lista_nodos_actual)
	    grilla_peq_osc.extend(grilla_nodos)
	    n_t = len(grilla_nodos)
	    lista_nodos_sgte = grilla_nodos[n_t - 1]
	    dif = dif_media_listas_nodos(lista_nodos_actual, lista_nodos_sgte)
	    lista_nodos_actual = lista_nodos_sgte
	    cont = cont + 1
	    print cont, dif
	equilibrio = lista_nodos_sgte
	n_t = len(grilla_peq_osc)
	i = 0
	dist_eq = amp_med_p_o + 1
	while(dist_eq > amp_med_p_o):
		dist_eq = dif_media_listas_nodos(grilla_peq_osc[i], equilibrio)
		i = i+1
	grilla_peq_osc = grilla_peq_osc[i: n_t - 1]
	return grilla_peq_osc


def dif_media_listas_nodos(lista_nodos_1, lista_nodos_2):
	spines_1 = sp.extraer_spines([lista_nodos_1])[:, 0, :]
	spines_2 = sp.extraer_spines([lista_nodos_2])[:, 0, :]
	n_nodos = np.size(spines_1, 0)
	dif = np.zeros(n_nodos)
	for i in range(n_nodos):
		dif_vector = np.abs(spines_2[i, :] - spines_1[i, :])
		dif[i] = sp.norma(dif_vector)
	dif_media = np.mean(dif)
	return dif_media


def H_sinusoidal(x, y, t, A=2, k=np.array([10, 0, 0]), w = 10, fase=0, h_unitario=np.array([0, 0, 1])):
	'''
	'''
	h_unitario = sp.normalizar(h_unitario)
	h = A * np.sin(k[0] * x + k[1] * y - w * t - fase)
	modulo_k = sp.norma(k)
	
	return h_unitario * h


def generar_arreglo_H(posiciones, t_array, A=2, k=np.array([10, 0, 0]), w = 10, fase=0, h_unitario=np.array([0, 0, 1])):
	'''
	'''
	n_nodos = np.size(posiciones, 0)
	n_t = np.size(t_array, 0)
	H_ext = np.zeros((n_t, n_nodos, 3))
	for i in range(n_t):
		for j in range(n_nodos):
			x = posiciones[j, 0]
			y = posiciones[j, 1]
			t = t_array[i]
			H_ext[i, j, :] = H_sinusoidal(x, y, t,  A=A, k=k, w = w, fase=fase, h_unitario=h_unitario)
	return H_ext



def extraer_fila(spines, nx, ny, j): 
	fila_spines = np.zeros((nx, 3), dtype=complex)
	indices=[]
	for i in range(nx):
		indice = indice_grilla_rectangular(nx, i, j)
		indices.append(indice)
		fila_spines[indice, :] = spines[indice, :]
	return fila_spines, indices



def extraer_columna(spines, nx, ny, i):
	columna_spines = np.zeros((ny, 3), dtype=complex)
	indices = []
	for j in range(ny):
		indice = indice_grilla_rectangular(nx, i, j)
		indices.append(indice)
		columna_spines[j, :] = spines[indice, :]
	return columna_spines, indices


def reemplazar_elementos(arreglo, elementos_reemplazo, indices):
	n = len(indices)
	for i in range(n):
		arreglo[indices[i]] = elementos_reemplazo[i]
	return arreglo

'''
def ft_gr_instante(spines, posiciones, nx, ny, d):
	ft_spines = np.zeros_like(spines, dtype=complex)
	ft_spines[:, :] = spines[:, :] 
	momentos = np.zeros_like(posiciones)
		
	for j in range(ny):
		# Transformar fila j.
		fila_spines, indices = extraer_fila(ft_spines, nx, ny, j)
		ft_fila_spines = np.zeros_like(fila_spines, dtype=complex)
		fila_kx = 2 * np.pi * fftfreq(nx, d)
		for k in range(3):
			ft_fila_spines[:, k] = fft(fila_spines[:, k]) / nx
		ft_spines = reemplazar_elementos(ft_spines, ft_fila_spines, indices)
		momentos[:, 0] = reemplazar_elementos(momentos[:, 0], fila_kx, indices)
		#print 'y=' + str(j), indices
		#print momentos
	
	for i in range(nx):
		# Transformar columna i
		columna_spines, indices = extraer_columna(ft_spines, nx, ny, i)
		ft_columna_spines = np.zeros_like(columna_spines, dtype=complex)
		columna_ky = 2 * np.pi * fftfreq(ny, d)
		for k in range(3):
			ft_columna_spines[:, k] = fft(columna_spines[:, k]) / ny
		ft_spines = reemplazar_elementos(ft_spines, ft_columna_spines, indices)
		momentos[:, 1] = reemplazar_elementos(momentos[:, 1], columna_ky, indices)
		#print 'x='+str(i), indices
		#print momentos

	#momentos[:, 1] = posiciones[:, 1]
	return np.abs(ft_spines), momentos


'''


def grilla_spines_instante(spines, nx, ny):
	grilla_spines = np.zeros((nx, ny, 3))
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, i, j)
			grilla_spines[j, i, :] = spines[indice_nodo, :]
	return grilla_spines


def ft_gr_instante(spines, nx, ny, d):
	grilla_spines = grilla_spines_instante(spines, nx, ny)
	ft_spines = np.zeros_like(grilla_spines)
	for i in range(3):
		modulo_ft = np.abs(fft2(grilla_spines[:, :, i]))
		ft_spines[:, :, i] = fftshift(modulo_ft)
	if not(is_odd(nx)):
		ft_spines = ft_spines[1:,:,:]
	if not(is_odd(ny)):
		ft_spines = ft_spines[:, 1:, :]
	return ft_spines

def is_odd(num):
    return num & 0x1


def onda_spin_sinusoidal(x, y, k=np.array([1, 1, 0])):
	arg = k[0] * x + k[1] * y
	sx = np.cos(arg)
	sy = np.sin(arg)
	sz = 0
	return np.array([sx, sy, sz])


def grilla_onda_spin(lista_nodos, posiciones, nx, ny, k):
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, i, j)
			[x, y] = posiciones[indice_nodo, :]
			s = onda_spin_sinusoidal(x, y, k=k)
			lista_nodos[indice_nodo].spin = s
	return lista_nodos



 


def recuperar_grilla_momentos(nx, ny, d=1.0):
	kx = 2 * np.pi * fftfreq(nx, d)
	ky = 2 * np.pi * fftfreq(ny, d)
	kx = fftshift(kx)
	ky = fftshift(ky)
	grilla_momentos = np.zeros((nx, ny, 2))
	for i in range(nx):
		for j in range(ny):
			grilla_momentos[i, j, :] = [kx[i], ky[j]]
	if not(is_odd(nx)):
		grilla_momentos = grilla_momentos[1:,:,:]
	if not(is_odd(ny)):
		grilla_momentos = grilla_momentos[:, 1:, :]
	return grilla_momentos



'''
def grilla_fourier(ft_spines, momentos, nx, ny, d=1.0):
	kx_array, ky_array = recuperar_frec_ordenadas(nx, ny, d=d)
	grilla_ft_spines = np.zeros((nx, ny, 3))
	for indice_nodo in range(nx * ny):
		[kx, ky] = momentos[indice_nodo, :]
		i = np.argmin(np.abs(kx_array - kx))
		j = np.argmin(np.abs(ky_array - ky))
		grilla_ft_spines[i, j, :] = ft_spines[indice_nodo, :] 
	return grilla_ft_spines

'''