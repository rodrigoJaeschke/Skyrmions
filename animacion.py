
# coding: utf-8
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import funciones_spin as sp
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation


def anima(indice_t, spines, posiciones):
    grilla_en_t = spines[:, indice_t, :]
    n_nodos = np.size(grilla_en_t, 0)
    
    n_lineas = 10
    sep_lineas = 2
    x_min, x_max = np.min(posiciones[:, 0]) -1.1, np.max(posiciones[:, 0]) +1.1
    y_min, y_max = np.min(posiciones[:, 1]) -1.1, np.max(posiciones[:, 1]) +1.1
    altura_plot = (x_max - x_min) / 2
    z_min, z_max = -altura_plot, altura_plot

    x = np.linspace(x_min, x_max, n_lineas)
    y = np.linspace(y_min, y_max, n_lineas)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros((n_lineas, n_lineas))
       
    ax = fig_anima.add_subplot(111, projection='3d')
    linea, = ax.plot([], [], [], label = 'Nodos', lw = 0.5)
    
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])

    #ax = Axes3D(fig)
    ax.azim = 80
    ax.elev = 60
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.set_zticklabels('')
    ax.plot_wireframe(X, Y, Z, rstride=sep_lineas, cstride=sep_lineas)
    for i in range(n_nodos):
        x0 = posiciones[i, 0]
        y0 = posiciones[i, 1]
        z0 = 0
        sx, sy, sz = grilla_en_t[i, 0], grilla_en_t[i, 1], grilla_en_t[i, 2]
        xf = x0 + sx
        yf = y0 + sy
        zf = z0 + sz
        color = sp.seleccionar_color(sx, sy, sz)
        arrow = sp.Arrow3D([x0, xf], [y0, yf], [z0, zf], mutation_scale=10, lw=2, arrowstyle="-|>", color=color)
        ax.add_artist(arrow)
        ax.plot([x0], [y0], [z0], color=color, markersize=3)
    return linea,


spines, posiciones = sp.cargar_spines()


fig_anima = plt.figure()
cuadros = np.size(spines, 1)
anim = animation.FuncAnimation(fig_anima, anima, frames = cuadros, fargs = (spines, posiciones), interval = 1, blit = True)


anim.save('img_animaciones/grillaboni.mp4', fps = 6)
