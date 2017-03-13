# coding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import funciones_spin as sp
import grilla_rectangular as gr
from scipy.fftpack import fft, fftfreq
from matplotlib.colors import LogNorm





prim_iter  = True
plotear = True

cmap = False
alternar = False
calcular_w_n = False
calcular_E = True

n = 60
nx = 35
ny = 20
[i_min, i_max, j_min, j_max] = [0, 35, 0, 20]



pi = np.pi


D = 1.0
j = -1.0
h = 0.36
A = 0.5
sigma= 11.07


'''
D = 0.6
j = 1.0
h = 0 #0.36
A = 0.1 #0.5
sigma=9.02
'''

print 'Obteniendo condicion inicial'
if not(prim_iter):
	spines, posiciones = sp.cargar_spines()
	n_t = np.size(spines, 1)
	condicion_inicial = spines[:, n_t-1, :]
	lista_nodos0, posiciones = gr.iniciar_grilla_rectangular(nx, ny, j, D, cb_periodica= False, modo='peq_osc', prim_iter=prim_iter, cond_inicial=condicion_inicial, interfacial=True)
	energias_antiferro = sp.energia_lista_nodos(lista_nodos0, h=h, A=A)

else:
	lista_nodos0, posiciones = gr.iniciar_grilla_rectangular(nx, ny, j, D, cb_periodica=False, modo='peq_osc', interfacial=True)
	energias_antiferro = sp.energia_lista_nodos(lista_nodos0, h=h, A=A)
	centro = np.array([nx, ny]) / 2
	#lista_nodos0 = gr.grilla_skyrmion(lista_nodos0, posiciones, nx, ny, sigma, gamma=0,  centro=centro)
	#lista_nodos0 = gr.alternar_lista_nodos(lista_nodos0, nx, ny)
#lista_nodos0 = gr.grilla_onda_spin(lista_nodos0, posiciones, n, n, np.array([2, 2, 0]))
#lista_nodos0 = gr.grilla_onda_triangular(lista_nodos0, posiciones, n, n, 1.0, phi=0)
#lista_nodos0 = gr.alternar_lista_nodos(lista_nodos0, n, n, epsilon=0.01)



'''
indice_x = 2
indice_y = 1
nodo_dmi = gr.indice_grilla_rectangular(nx, ny, indice_x, indice_y)
lista_nodos0 = gr.invertir_spin(lista_nodos0, nx, ny, indice_x, indice_y)
'''


'''
k_mod = 2 * np.pi / n
#k_unit = np.array([1, 1, 0)
#k_unit = sp.normalizar(k_unit)

k = k_mod * k_unit

w = 10
A = 2
'''

print 'Simulando avance temporal'
h_unitario = np.array([-1, 1, 0])
step = 0.1
tmax = 0.2
n_steps = int(tmax/step)
t_array = np.linspace(0, tmax, n_steps)
#H_ext = gr.generar_arreglo_H(posiciones, t_array, k=np.array([0, 0, 0]), w=0, A = 0, h_unitario=h_unitario)
H_ext = gr. generar_arreglo_H_cte(posiciones, t_array, np.array([0, 0, h]))
grilla_nodos = sp.simular_spines(t_array, lista_nodos0, H_ext=H_ext, A=A)
spines = sp.extraer_spines(grilla_nodos)

'''
print 'Guardando grilla de spines'
sp.guardar_spines(spines, posiciones)
print 'Grilla de spines guardada'
'''

cuadro = 0
eje_spin = 0

cuadro2 = np.size(t_array, 0) - 1



'''
for i in range(len(grilla_nodos)):
	for k in range(len(grilla_nodos[i])):
		print k
		print grilla_nodos[i][k].spin
		print grilla_nodos[i][k].vecinos
'''


if plotear:
	spines_en_t = spines[:, cuadro, :]
	ft_spines = gr.ft_gr_instante(spines_en_t, nx, ny, 1.0)
	momentos = gr.recuperar_grilla_momentos(n, n)
	limites_imshow = [0, 2 * pi, 0, 2 * pi]
	ft_imagen = ft_spines[:, ::-1, eje_spin].transpose()
	fig1 = plt.figure(1)
	plt.imshow(ft_imagen, extent=limites_imshow, cmap='viridis')
	plt.plot([0, 2*pi], [pi, pi], linewidth=2, color='black')
	plt.plot([pi, pi], [0, 2*pi], linewidth=2, color='black')
	plt.xlim(0, 2*pi)
	plt.ylim(0, 2*pi)
	plt.title('Transformada de fourier de la grilla de spines.')
	plt.colorbar()
	plt.grid()
	plt.xlabel('$k_x$', fontsize=15)
	plt.ylabel('$k_y$', fontsize=15)

	if cmap:
		fig2 = gr.plotear_instante_grilla_cmap(spines, cuadro, posiciones, (nx, ny), 2, show=False)
		fig3 = gr.plotear_instante_grilla_cmap(spines, cuadro2, posiciones, (nx, ny), 3, show=False)
	else:
		mask = gr.mascara_rectangular(nx, ny, i_min, i_max, j_min, j_max)
		fig2 = plt.figure(2)
		sp.plotear_instante_grilla(spines, cuadro, posiciones, fig2, mask=mask)
		plt.title('primer_cuadro')


		fig3 = plt.figure(3)
		sp.plotear_instante_grilla(spines, cuadro2, posiciones, fig3, mask=mask)
		plt.title('ultimo cuadro')









path = 'Graficos/'#'/Documentos/Universidad/6to_semestre/Skyrmions/Skyrmions/Graficos/'
terminacion= '.png'
nombre_fig1 = 'j'+str(j)+'D'+str(D)+'H'+str(h)+'fourier'+terminacion
nombre_fig2 = 'j'+str(j)+'D'+str(D)+'H'+str(h)+'cuadro'+str(cuadro)+terminacion
nombre_fig3 = 'j'+str(j)+'D'+str(D)+'H'+str(h)+'cuadro'+str(cuadro2)+terminacion



'''
fig1.savefig(path+nombre_fig1)
fig2.savefig(path+nombre_fig2)
fig3.savefig(path+nombre_fig3)
'''


if calcular_w_n:
	w_n = gr.wending_number_rectangular(spines[:, cuadro2, :], nx, ny)
	print 'Wending Number = ', w_n


if calcular_E:
	n_t = len(grilla_nodos)
	
	energias = np.zeros((n_t, 5))
	for i in range(n_t):
		energias[i, :] = sp.energia_lista_nodos(grilla_nodos[i], h=h, A=A)

	plt.figure(4)
	plt.plot(t_array, energias[:, 0], '-')
	plt.title('Energia en funcion del tiempo')
	plt.xlabel('tiempo')
	plt.ylabel('Energia total')

	print energias[0, 0]
	print energias[n_t-1, 0]
	

	if False:
		energia_0 = sp.energia_lista_nodos(grilla_nodos[0])
		print 'Energia inicial = ', energia_0
		print ' '
		print ' '
		print ' '
		energia_final = sp.energia_lista_nodos(grilla_nodos[n_t - 1])
		print 'Energia final = ', energia_final
	

plt.show()