# coding: utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import funciones_spin as sp
import grilla_rectangular as gr
from scipy.fftpack import fft, fftfreq
from matplotlib.colors import LogNorm

pi = np.pi




n = 20
lista_nodos0, posiciones = gr.iniciar_grilla_rectangular(n, n, 0.7, modo='peq_osc')
lista_nodos0 = gr.grilla_onda_spin(lista_nodos0, posiciones, n, n, np.array([-1, 1, 0]))
k_mod = 2 * np.pi / n
k_unit = np.array([1, 1, 0])
k_unit = sp.normalizar(k_unit)

k = k_mod * k_unit
h_unitario = np.array([-1, 1, 0])
w = 10
A = 2


step = 0.1
tmax = 10
n_steps = int(tmax/step)
t_array = np.linspace(0, tmax, n_steps)
H_ext = gr.generar_arreglo_H(posiciones, t_array, k=k, w=w, A = A, h_unitario=h_unitario)
grilla_nodos = sp.simular_spines(t_array, lista_nodos0, H_ext=0)
spines = sp.extraer_spines(grilla_nodos)
sp.guardar_spines(spines, posiciones)



cuadro = 60
eje_spin = 1

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
spines_en_t = spines[:, cuadro, :]
ft_spines = gr.ft_gr_instante(spines_en_t, n, n, 1.0)
momentos = gr.recuperar_grilla_momentos(n, n)


nx = np.size(momentos, 0)
ny = np.size(momentos, 1)

for i in range(nx):
	for j in range(ny):
		ax.scatter(momentos[i, j, 0], momentos[i, j, 1], ft_spines[i, j, eje_spin])
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('($S_x$)')




limites_imshow = [-pi, pi, pi, -pi]
plt.figure(2)
plt.imshow(ft_spines[:, :, eje_spin], extent=limites_imshow, cmap='jet')
plt.title('Transformada de fourier de la grilla de spines.')
plt.colorbar()
plt.grid()
plt.ylabel('$k_x$')
plt.xlabel('$k_y$')




fig4 = plt.figure(4)
ax = fig4.add_subplot(111, projection='3d')

lista_ft_spines = np.zeros((nx * ny, 3))
lista_momentos = np.zeros((nx * ny, 2))
for i in range(nx):
	for j in range(ny):
		indice_nodo = gr.indice_grilla_rectangular(nx, i, j)
		lista_ft_spines[indice_nodo, :] = ft_spines[i, j, :]
		lista_momentos[indice_nodo, :] = momentos[i, j, :]

ax.plot_surface(lista_momentos[:, 0], lista_momentos[:, 1], lista_ft_spines[:, eje_spin], cmap='hot')



fig3 = plt.figure(3)
sp.plotear_instante_grilla(spines, cuadro, posiciones, fig3)

'''
for i in range(nx):
	for j in range(ny):
		print ft_spines[i, j, 0]
'''

plt.show()

