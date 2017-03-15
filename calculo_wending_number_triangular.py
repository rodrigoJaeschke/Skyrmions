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


n = 60
nx = n
ny = n

spines, posiciones = sp.cargar_spines()
n_t = np.size(spines, 1)
w_n = np.zeros(n_t)
for i in range(n_t):
	spines_en_t = spines[:, i, :]
	w_n[i] = gt.wending_number_triangular(spines_en_t, nx, ny)

print w_n[0]

plt.plot(range(n_t), w_n)
plt.show()