# coding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import funciones_spin as sp
import grilla_rectangular as gr
from scipy.fftpack import fft, fftfreq
from matplotlib.colors import LogNorm
import grilla_triangular as gt



prim_iter = False
plotear = True
calcular_E = True

[i_min, i_max, j_min, j_max] = [31, 39, 33, 42]



n = 60
nx = n
ny = n


'''
D = 1.0
j = 0
h = 0
A = 0
'''


D = 1
j = -1
h = 0.36
A = 0.5


print 'Obteniendo condicion inicial'
if not(prim_iter):
	spines, posiciones = sp.cargar_spines()
	n_t = np.size(spines, 1)
	condicion_inicial = spines[:, n_t - 1, :]
	lista_nodos0, posiciones = gt.iniciar_grilla_triangular(nx, ny, j, D, cb_periodica=True, prim_iter=prim_iter, cond_inicial = condicion_inicial)	

else:
	lista_nodos0, posiciones = gt.iniciar_grilla_triangular(nx, ny, j, D, cb_periodica=True, prim_iter=prim_iter)


print 'Simulando avance temporal'

step = 0.01
tmax = 0.02
n_steps = int(tmax/step)
t_array = np.linspace(0, tmax, n_steps)

H_ext = gr.generar_arreglo_H_cte(posiciones, t_array, np.array([0, 0, h]))
grilla_nodos = sp.simular_spines(t_array, lista_nodos0, H_ext=H_ext, A=A)
spines = sp.extraer_spines(grilla_nodos)

print 'Guardando grilla de spines'
sp.guardar_spines(spines, posiciones)
print 'Grilla de spines guardada'



cuadro = 0
eje_spin = 0
cuadro2 = np.size(t_array, 0) - 1




#fig2 = gr.plotear_instante_grilla_cmap(spines, cuadro, posiciones, (nx, ny), 2, show=False)
#fig3 = gr.plotear_instante_grilla_cmap(spines, cuadro2, posiciones, (nx, ny), 3, show=False)


if plotear:
	mask = gr.mascara_rectangular(nx, ny, i_min, i_max, j_min, j_max)
	fig2 = plt.figure(2)
	sp.plotear_instante_grilla(spines, cuadro, posiciones, fig2, mask=mask)
	plt.title('primer_cuadro')

	fig3 = plt.figure(3)
	sp.plotear_instante_grilla(spines, cuadro2, posiciones, fig3, mask=mask)
	plt.title('Acercamiento de un Skyrmion ferromagnetico')




'''
n_t = len(grilla_nodos)
energias = np.zeros(n_t)
for i in range(n_t):
	energias[i] = sp.energia_lista_nodos(grilla_nodos[i], A)


plt.figure(4)
plt.plot(t_array, energias, '-')
plt.title('Energia en funcion del tiempo')
plt.xlabel('tiempo')
plt.ylabel('Energia total')

'''

if calcular_E:
	n_t = len(grilla_nodos)
	
	energias = np.zeros((n_t, 5))
	for i in range(n_t):
		energias[i, :] = sp.energia_lista_nodos(grilla_nodos[i], h, A)

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