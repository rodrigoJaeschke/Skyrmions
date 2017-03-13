# coding: utf-8
from __future__ import division
import numpy as np
import funciones_spin as sp
from scipy.fftpack import fft, fftfreq, fft2, fftshift
import random as rand
import matplotlib.pyplot as plt

def iniciar_grilla_rectangular(nx, ny, j, D, modo='aleatorio', cb_periodica=True, cond_inicial=0, prim_iter=True, interfacial=False):
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
	listas_vecinos = lista_vecinos_rectangular(nx, ny, j, D, cb_periodica, interfacial)
	if prim_iter:
		if (modo == 'peq_osc'):
			print 'Condicion inicial ferro/antiferro (dependiendo del signo de j)'
			condicion_inicial= equilibrio_rectangular(nx, ny, j)
		else:
			print 'Condicion inicial aleatoria.'
			condicion_inicial = condicion_inicial_aleatoria(nx * ny)
	else:
		print 'recuperando CI anterior'
		condicion_inicial = cond_inicial
	lista_nodos = []
	for l in range(ny):
		for k in range(nx):
			i = indice_grilla_rectangular(nx, ny, k, l)
			posiciones[i, :] = [X[l, k], Y[l, k]]
			nodo_nuevo = sp.Nodo_spin(condicion_inicial[i], H_ext, listas_vecinos[i])
			lista_nodos.append(nodo_nuevo)
	return lista_nodos, posiciones


def perturbacion(epsilon):
	ds_x = (rand.random() - 0.5) * epsilon 
	ds_y = (rand.random() - 0.5) * epsilon
	ds_z = (rand.random() - 0.5) * epsilon
	ds = np.array([ds_x, ds_y, ds_z])
	return ds


def equilibrio_rectangular(nx, ny, j, s_eq=np.array([0, 0, 1]), epsilon=0.01):
	cond_inicial = np.zeros((nx * ny, 3))
	s_eq = sp.normalizar(s_eq)
	for k in range(nx):
		for l in range(ny):
			i = indice_grilla_rectangular(nx, ny, k, l)
			ds = perturbacion(epsilon)
			if (j > 0):
				s = sp.normalizar(s_eq * (-1) ** (k + l) + ds)
			else:
				s = sp.normalizar(s_eq + ds)
			cond_inicial[i] = s
	return cond_inicial

def coord_grilla_rectangular(i, nx):
	indice_y = int(i / nx)
	indice_x = i - indice_y * nx
	return indice_x, indice_y

def D_ij_rectangular(D, indice_nodo, indice_vecino, nx, ny, interfacial):
	x_nodo, y_nodo = coord_grilla_rectangular(indice_nodo, nx)
	x_vec, y_vec = coord_grilla_rectangular(indice_vecino, nx)
	Dx = x_vec - x_nodo
	Dy = y_vec - y_nodo	
	D_vector = D * sp.normalizar(np.array([Dx, Dy, 0]))
	borde_izq = (x_nodo == 0) and (x_vec == nx - 1)
	borde_der = (x_nodo == nx - 1) and (x_vec == 0)
	borde_sup = (y_nodo == 0) and (y_vec == ny - 1)
	borde_inf = (y_nodo == ny - 1) and (y_vec == 0)
	cruza_borde = borde_izq or borde_der or borde_inf or borde_sup
	if cruza_borde:
		D_vector = - D_vector
	if interfacial:
		D_vector[0:2] = rotar_2d(D_vector[0:2], np.pi / 2)

	return D_vector


def lista_vecinos_rectangular(nx, ny, j_intercambio, D, cb_periodica, interfacial):
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
	for j in range(ny):
	    for i in range(nx):
	    	indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
	        indices = []
	        if not(j == 0):
	            indices.append(indice_grilla_rectangular(nx, ny, i, j - 1))
	        elif cb_periodica:
	            indices.append(indice_grilla_rectangular(nx, ny, i, -1))
	            
	        if not(i == nx - 1):
	            indices.append(indice_grilla_rectangular(nx, ny, i + 1, j))
	        elif cb_periodica:
	            indices.append(indice_grilla_rectangular(nx, ny, 0, j))
	            
	        if not(j == ny - 1):
	            indices.append(indice_grilla_rectangular(nx, ny, i, j + 1))
	        elif cb_periodica:
	            indices.append(indice_grilla_rectangular(nx, ny, i, 0))
	            
	        if not(i == 0):
	            indices.append(indice_grilla_rectangular(nx, ny, i - 1, j))
	        elif cb_periodica:
	            indices.append(indice_grilla_rectangular(nx, ny, -1, j))

	        #print 'indice Nodo: ', indice_nodo
	        #print 'indices_vecinos: ', indices
	            
	        vecinos = np.zeros((len(indices), 5))
	        for k in range(len(indices)):
	        	#print 'indice vecino: ', indices[k]
	        	vecinos[k, 0:2] = [indices[k], j_intercambio]
	        	D_vector = D_ij_rectangular(D, indice_nodo, indices[k], nx, ny, interfacial)
	        	#print 'D = ' , D_vector
	        	vecinos[k, 2:5] = D_vector
	        listas_vecinos.append(vecinos)
	return listas_vecinos


def indice_grilla_rectangular(nx, ny, i, j):
	'''
	Recibe un taman~o nx (horizontal) de una grilla
	rectangular de nodos, asi como los indices
	(i, j) que definen las coordenadas de un determinado
	Nodo_spin en la grlla.

	La función retorna el índice del correspondiente nodo
	(i, j), con el que se enumerarian los nodos en un arreglo
	unidimensional.
	'''
	if i == -1:
		i = nx - 1
	if j == -1:
		j = ny - 1
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

def generar_arreglo_H_cte(posiciones, t_array, h):
	n_nodos = np.size(posiciones, 0)
	n_t = np.size(t_array, 0)
	H_ext = np.zeros((n_t, n_nodos, 3))
	for i in range(n_t):
		for j in range(n_nodos):
			H_ext[i, j, :] = h
	return H_ext



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






def grilla_spines_instante(spines, nx, ny, dim=3):
	grilla_spines = np.zeros((nx, ny, dim))
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
			grilla_spines[i, j, :] = spines[indice_nodo, :]
	return grilla_spines


def ft_gr_instante(spines, nx, ny, d):
	grilla_spines = grilla_spines_instante(spines, nx, ny)
	ft_spines = np.zeros_like(grilla_spines)
	for i in range(3):
		modulo_ft = np.abs(fft2(grilla_spines[:, :, i]))
		ft_spines[:, :, i] = modulo_ft #fftshift(modulo_ft)



	ft_spines_completado = np.zeros((nx + 1, ny + 1, 3))
	ft_spines_completado[0:nx, 0:ny, :] = ft_spines[:, :, :]
	ft_spines_completado[nx, :, :] = ft_spines_completado[0, :, :]
	ft_spines_completado[:, ny, :] = ft_spines_completado[:, 0, :]
	return ft_spines_completado

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
			indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
			[x, y] = posiciones[indice_nodo, :]
			s = onda_spin_sinusoidal(x, y, k=k)
			lista_nodos[indice_nodo].spin = s
	return lista_nodos



 


def recuperar_grilla_momentos(nx, ny, d=1.0):
	kx = 2 * np.pi * fftfreq(nx, d)
	ky = 2 * np.pi * fftfreq(ny, d)
	#kx = fftshift(kx)
	#ky = fftshift(ky)
	grilla_momentos = np.zeros((nx, ny, 2))
	for i in range(nx):
		for j in range(ny):
			grilla_momentos[i, j, :] = [kx[i], ky[j]]
	'''
	if not(is_odd(nx)):
		grilla_momentos = grilla_momentos[1:,:,:]
	if not(is_odd(ny)):
		grilla_momentos = grilla_momentos[:, 1:, :]
	'''
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

def alternar_lista_nodos(lista_nodos, nx, ny, epsilon=0):
	for i in range(nx):
		for j in range(ny):
			indice = indice_grilla_rectangular(nx, ny, i, j)
			nodo = lista_nodos[indice]
			ds = perturbacion(epsilon)
			nodo.spin = sp.normalizar((-1)**(i+j) * nodo.spin + ds)
	return lista_nodos


def onda_triangular(x, y, modulo_k, phi):
	a = 0.5
	b = np.sqrt(3) / 2
	k1 = modulo_k * np.array([1, 0])
	k2 = modulo_k * np.array([-a, b])
	k3 = modulo_k * np.array([-a, -b])

	k1 = rotar_2d(k1, phi)
	k2 = rotar_2d(k2, phi)
	k3 = rotar_2d(k3, phi)

	sx1 = np.cos(k1[0] * x + k1[1] * y)
	sx2 = np.cos(k2[0] * x + k2[1] * y)
	sx3 = np.cos(k3[0] * x + k3[1] * y)

	sy1 = np.sin(k1[0] * x + k1[1] * y)
	sy2 = np.sin(k2[0] * x + k2[1] * y)
	sy3 = np.sin(k3[0] * x + k3[1] * y)

	sx = sx1 + sx2 + sx3
	sy = sy1 + sy2 + sy3
	sz = 0
	s = sp.normalizar(np.array([sx, sy, sz]))
	return s

def grilla_onda_triangular(lista_nodos, posiciones, nx, ny, modulo_k, phi=0):
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
			[x, y] = posiciones[indice_nodo, :]
			s = onda_triangular(x, y, modulo_k, phi)
			lista_nodos[indice_nodo].spin = s
	return lista_nodos


def rotar_2d(vector_2d, phi):
	mat_rot =np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])
	vec_rotado = np.matmul(mat_rot, vector_2d)
	return vec_rotado

def cart_to_pol(x, y):
	r = np.sqrt(x ** 2 + y ** 2)
	phi = np.arctan2(y, x)
	return r, phi


def funcion_skyrmion(x, y, x0, y0, sigma, m, gamma):
	x = x - x0
	y = y - y0
	r, phi = cart_to_pol(x, y)
	theta_spin = np.pi * np.exp(-r ** 2 / (2 * sigma ** 2))
	phi_spin = m * phi + gamma
	spin = sp.esf_to_cart(1, theta_spin, phi_spin)
	return spin

def grilla_skyrmion(lista_nodos, posiciones, nx, ny, sigma, m=1, gamma=0, centro=np.array([0, 0])):
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
			[x, y] = posiciones[indice_nodo, :]
			s = funcion_skyrmion(x, y, centro[0], centro[1], sigma, m, gamma)
			lista_nodos[indice_nodo].spin = s
	return lista_nodos


def invertir_spin(lista_nodos, nx, ny, i, j):
	indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
	nodo = lista_nodos[indice_nodo]
	nodo.spin = -1 * nodo.spin
	return lista_nodos



def plotear_grilla_con_dmi(grilla_nodos, posiciones, nodo_dmi, indice_vec, fig, indice_t=0):
	spines = sp.extraer_spines(grilla_nodos)
	grilla_en_t = spines[:, indice_t, :]
	n_nodos = np.size(grilla_en_t, 0)
	ax = fig.add_subplot(111, projection='3d')
	x_min, x_max = np.min(posiciones[:, 0]) -1.1, np.max(posiciones[:, 0]) +1.1
	y_min, y_max = np.min(posiciones[:, 1]) -1.1, np.max(posiciones[:, 1]) +1.1
	altura_plot = (x_max - x_min) / 2
	z_min, z_max = -altura_plot, altura_plot
	ax.set_xlim([x_min, x_max])
	ax.set_ylim([y_min, y_max])
	ax.set_zlim([z_min, z_max])
	n_lineas = 10
	sep_lineas = 2
	x = np.linspace(x_min, x_max, n_lineas)
	y = np.linspace(y_min, y_max, n_lineas)
	X, Y = np.meshgrid(x, y)
	Z = np.zeros((n_lineas, n_lineas))
	ax.plot_wireframe(X, Y, Z, rstride=sep_lineas, cstride=sep_lineas)
	ax.azim = 80
	ax.elev = 60
	for i in range(n_nodos):
		x0 = posiciones[i, 0]
		y0 = posiciones[i, 1]
		z0 = 0
		sx, sy, sz = grilla_en_t[i, 0], grilla_en_t[i, 1], grilla_en_t[i, 2]
		xf = x0 + sx
		yf = y0 + sy
		zf = z0 + sz
		color = sp.seleccionar_color(sx, sy, sz)
		arrow = sp.Arrow3D([x0, xf], [y0, yf], [z0, zf], mutation_scale=10, lw=3, arrowstyle="-|>", color=color)
		ax.add_artist(arrow)
		ax.plot([x0], [y0], [z0], '.', color=(0, 0, 0), markersize=5)
		if i==nodo_dmi:
			vecinos_i = grilla_nodos[indice_t][i].vecinos
			n_vec = np.size(vecinos_i, 0)
			for j in range(n_vec):
				if j==indice_vec:
					x0 = posiciones[i, 0]
					y0 = posiciones[i, 1]
					Dx = vecinos_i[j, 2]
					Dy = vecinos_i[j, 3]
					xf = x0 + Dx
					yf = y0 + Dy
					arrow = sp.Arrow3D([x0, xf], [y0, yf], [0, 0], mutation_scale=10, lw=3, arrowstyle="-|>", color='black')
					ax.add_artist(arrow)
					print 'vecino'+str(j)+' del nodo '+str(i)+', DMI= ', [Dx, Dy, 0]

	ax.set_xlabel('$X$')
	ax.set_ylabel('$Y$')
	ax.set_zlabel('$Z$')
	return ax


def posiciones_3x3():
	posiciones = np.zeros((9, 2))
	for i in range(3):
		for j in range(3):
			indice = indice_grilla_rectangular(3, 3, i, j)
			posiciones[indice, :] = [i, j]
	return posiciones



def plotear_instante_grilla_cmap(spines, indice_t, posiciones, shape, num_fig, show=False):
	(nx, ny) = shape
	spines_en_t = spines[:, indice_t, :]
	grilla_spines = grilla_spines_instante(spines_en_t, nx, ny)
	grilla_posiciones = grilla_spines_instante(posiciones, nx, ny, dim=2)
	spines_z = grilla_spines[:, :, 2]
	spines_x = grilla_spines[:, :, 0]
	spines_y = grilla_spines[:, :, 1]
	X = grilla_posiciones[:, :, 0]
	Y = grilla_posiciones[:, :, 1]
	fig = plt.figure(num_fig)	
	imagen = spines_z[:, ::-1].transpose()
	
	x0, x_max = np.min(posiciones[:, 0]), np.max(posiciones[:, 0])
	y0, y_max = np.min(posiciones[:, 1]), np.max(posiciones[:, 1])
	limites = [x0, x_max, y0, y_max]

	title = 'Grilla de spines en el cuadro: ' + str(indice_t)
	X =np.append(X, [0])
	Y =np.append(Y, [0])
	spines_x =np.append(spines_x, [1.0])
	spines_y =np.append(spines_y, [0])
	plt.imshow(imagen, extent=limites, cmap='jet', vmin=-1, vmax=1, interpolation='nearest')
	plt.colorbar()
	Q = plt.quiver(X, Y, spines_x, spines_y, pivot='tail')
	qk = plt.quiverkey(Q, 0.1, 1.1, 1, r'1', labelpos='W',
                  fontproperties={'weight': 'bold'})

	plt.xlabel('$X$')
	plt.ylabel('$Y$')
	plt.title(title)
	plt.grid()
	return fig
	



def mascara_rectangular(nx, ny, i_min, i_max, j_min, j_max):
		mask = []
		for j in range(ny):
			for i in range(nx):
				indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
				if (i >= i_min) and (i < i_max) and (j >= j_min) and (j < j_max):
					mask.append(indice_nodo)
		return mask



def wending_number_rectangular(spines_en_t, nx, ny):
	angulo_solido = 0
	for j in range(ny):
		for i in range(nx):
			area_triangulos = w_n_nodo_rectangular(spines_en_t, nx, ny, i, j)
			angulo_solido = angulo_solido + area_triangulos

	return angulo_solido / 4 / np.pi




def w_n_nodo_rectangular(spines_en_t, nx, ny, i, j):
	if not(i == nx - 1):
		indice_0 = indice_grilla_rectangular(nx, ny, i + 1, j)
	else:
	 	indice_0= indice_grilla_rectangular(nx, ny, 0, j)


	if not((i == nx - 1) or (j == ny - 1)):
		indice_1 = indice_grilla_rectangular(nx, ny, i + 1, j + 1)
	elif (i == nx - 1 and not(j == ny - 1)):
		indice_1 = indice_grilla_rectangular(nx, ny, 0, j + 1)
	elif not (i == nx - 1) and (j == ny - 1):
		indice_1 = indice_grilla_rectangular(nx, ny, i + 1, 0)
	else:
		indice_1 = indice_grilla_rectangular(nx, ny, 0, 0)

	if not(j == ny - 1):
		indice_2 = indice_grilla_rectangular(nx, ny, i, j + 1)
	else:
		indice_2 = indice_grilla_rectangular(nx, ny, i, 0)

	s_nodo = spines_en_t[indice_grilla_rectangular(nx, ny, i, j)]
	vertice_0 = spines_en_t[indice_0]
	vertice_1 = spines_en_t[indice_1]
	vertice_2 = spines_en_t[indice_2]
	area_triangulo_1 = sp.solid_angle_triangle(s_nodo, vertice_0, vertice_1)
	area_triangulo_2 = sp.solid_angle_triangle(s_nodo, vertice_1, vertice_2)

	return area_triangulo_1 + area_triangulo_2



'''
def reconstruir_grilla_nodos(spines, nx, ny, j, D):
	n
'''



def funcion_skyrmion2(x, y, x0, y0, R, W, m=1, gamma=0):
	x = x - x0
	y = y - y0
	r, phi = cart_to_pol(x, y)
	theta_spin = np.pi * np.arctan(np.exp(-(r-R)/W)) / np.arctan(np.exp(R/W))
	phi_spin = m * phi + gamma
	spin = sp.esf_to_cart(1, theta_spin, phi_spin)
	return spin

def grilla_skyrmion2(lista_nodos, posiciones, nx, ny, R, W, m=1, gamma=0, centro=np.array([0, 0])):
	for i in range(nx):
		for j in range(ny):
			indice_nodo = indice_grilla_rectangular(nx, ny, i, j)
			[x, y] = posiciones[indice_nodo, :]
			s = funcion_skyrmion2(x, y, centro[0], centro[1], sigma, m, gamma)
			lista_nodos[indice_nodo].spin = s
	return lista_nodos
