# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as lin
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation

class Nodo_spin:
    alpha = 0.5
    
    def __init__(self, s, H,  vecinos):
        self.campo_externo = H
        self.vecinos = vecinos
        self.spin = s
        

def set_alpha(lista_nodos, alpha):
    n_nodos = len(lista_nodos)
    for i in range(n_nodos):
        lista_nodos[i].alpha = alpha


def set_H_ext(lista_nodos, H_ext):
    n_nodos = len(lista_nodos)
    for i in range(n_nodos):
        lista_nodos[i].campo_externo = H_ext[i, :]
    return lista_nodos


def h_intercambio(lista_nodos, i):
    nodo = lista_nodos[i]
    vecinos = nodo.vecinos
    n_vec = np.size(vecinos, 0)
    h = np.zeros((3))
    for j in range(n_vec):
        indice_vec = int(vecinos[j, 0])
        s_vec = lista_nodos[indice_vec].spin
        j_vec = vecinos[j, 1]
        h = h - s_vec*j_vec
    return h


def resolver_LL(dt, H_eff, s_0):
    '''
    Recibe un salto de tiempo, el campo magnetico
    efectivo, y un vector de spin inicial s0.
    Retorna la siguiente iteracion del arreglo spin,
    precesandolo en torno al campo.
    '''
    I = np.zeros((3,3)) # matriz identidad
    I[0,0], I[1,1], I[2,2] = 1, 1, 1
    s = np.zeros((3))
    A = I - 0.5 * Lh(H_eff) * dt
    B = I + 0.5 * Lh(H_eff) * dt
    b = np.matmul(B, s_0)
    detA = det_dim3(A)
    for j in range(3):
        Aj = reemplazo_cramer(A, b, j)
        detAj = det_dim3(Aj)
        s[j] = detAj
    s = s / detA
    return s


def avanzar_nodo(dt, lista_nodos, i):
    '''
    Recibe el valor de la grilla de nodos,
    en un instante determinado, en una lista
    de nodos. Avanza el nodo de indice i.
    calculando el campo efectivo que siente
    producto del campo externo y las interacciones
    con sus vecinos.
    Retorna un nuevo nodo "avanzado" un dt en la Ec de mov.
    '''
    nodo = lista_nodos[i]
    spin = nodo.spin
    alpha = nodo.alpha
    h_ext = nodo.campo_externo
    h_int = h_intercambio(lista_nodos, i)
    H = h_ext + h_int
    H_eff =  H - alpha * cruz(spin, H)
    s = resolver_LL(dt, H_eff, spin)
    nuevo_nodo = Nodo_spin(s, h_ext, nodo.vecinos)
    return nuevo_nodo

def avanzar_grilla(dt, lista_nodos):
    '''
    Recibe una grilla con los nodos en un determinado
    instante, calcula los objetos nodo en la siguiente
    iteracion, y los retorna en una nueva lista de nodos.
    '''
    n = len(lista_nodos)
    nueva_lista_nodos = []
    for i in range(n):
        nuevo_nodo = avanzar_nodo(dt, lista_nodos, i)
        nueva_lista_nodos.append(nuevo_nodo)
    return nueva_lista_nodos

def simular_spines(t_array, lista_nodos0, H_ext = 0):
    '''
    Recibe un arreglo de tiempos, y una lista
    de nodos del instante inicial.
    Además recibe como parámetro opcional un arreglo,
    con los campos externos en cada nodo para cada tiempo,
    ordenados en formato:
    H_ext[indice_tiempo, indice_nodo, componente_H]
    (arreglo de 3 dimenciones)
    La funcion retorna una matriz de 2 dimensiones
    conteniendo las grillas de nodos
    en cada uno de los instantes de t_array.
    El formato de retorno es:
    [instante_tiempo, nodo_i]
    '''
    n_t = np.size(t_array,0)
    n_nodos = len(lista_nodos0)
    if type(H_ext) == int:
        H_ext = np.zeros((n_t, n_nodos, 3))
    set_H_ext(lista_nodos0, H_ext[0])
    lista_nodos = [lista_nodos0]
    for i in range(1, n_t):
        dt = t_array[i]- t_array[i-1]
        nueva_lista = avanzar_grilla(dt, lista_nodos[i-1])
        set_H_ext(nueva_lista, H_ext[i, :, :])
        lista_nodos.append(nueva_lista)
    return lista_nodos

def extraer_spines(lista_nodos):
    '''
    Recibe una matriz de nodos,en formato
    [Instane_tiempo, _nodo_i] y extrae los valores
    de los spines, en formato:
    [Nodo_i, Instante_tiempo, componentes_spin]
    como un arreglo de numpy.
    '''
    n_t = len(lista_nodos)
    n_nodos = len(lista_nodos[0])
    spines = np.zeros((n_nodos, n_t, 3))
    for i in range(n_t):
        for j in range(n_nodos):
            spines[j, i, :] = lista_nodos[i][j].spin
    return spines


    



def generar_latitudes(n):
    n_phi = 100
    theta = np.linspace(0, np.pi, n)
    phi = np.linspace(0, 2*np.pi, n_phi)
    curvas = np.zeros((n, n_phi, 3))
    for i in range(n):
        for j in range(n_phi):
            curvas[i, j, :] = esf_to_cart(1, theta[i], phi[j])
    return curvas

def generar_meridianos(n):
    n_theta = 100
    theta = np.linspace(0, 2*np.pi, n_theta)
    phi = np.linspace(0, 2*np.pi, n)
    curvas = np.zeros((n, n_theta, 3))
    for i in range(n):
        for j in range(n_theta):
            curvas[i, j, :] = esf_to_cart(1, theta[j], phi[i])
    return curvas

def plot_spin(s, t, figura):
    fig = plt.figure(figura)
    ax = fig.gca(projection='3d')
    ax.plot(s[:,0],s[:,1], zs=s[:,2], color='red')
    ax.plot([0],[0],[0],'.r')
    curvas_l = generar_latitudes(5)
    n_latitudes = 5
    n_meridianos = 8
    for i in range(n_latitudes):
        curva_i = curvas_l[i,:,:]
        x = curva_i[:,0]
        y = curva_i[:,1]
        z = curva_i[:,2]
        ax.plot(x, y, zs = z, marker='.', markersize=1, linewidth=0.1, color='blue')
    curvas_m = generar_meridianos(n_meridianos)
    for i in range(n_meridianos):
        curva_i = curvas_m[i,:,:]
        x = curva_i[:,0]
        y = curva_i[:,1]
        z = curva_i[:,2]
        ax.plot(x, y, zs = z, marker='.', markersize=1, linewidth = 0.1, color='blue')
    ax.set_xlabel('$X$', fontsize=15)
    ax.set_ylabel('$Y$', fontsize=15)
    ax.set_zlabel('$Z$', fontsize=15)
    plt.show()
    return 0

def punto(x, y):
    xy = x*y
    suma = 0
    for i in range(np.size(x)):
        suma = suma + xy[i]
    return suma    
    
    
def cruz(x,y):
    producto=np.zeros(3)
    producto[0] = x[1]*y[2] - x[2]*y[1]
    producto[1] = -(x[0]*y[2]-x[2]*y[0])
    producto[2] = x[0]*y[1] -x[1]*y[0]
    return producto

def norma(x):
    x2=x*x
    norma=0
    for i in range(0,np.size(x)):
        norma = norma+x2[i]
    norma = np.sqrt(norma)
    return norma

def normalizar(x):
    norma_x=norma(x)
    return x/norma_x


def Lh(h):
    L = np.zeros((3, 3, 3)) #tensor de rotacion, L[i,:,:] = L_i
    L[0, 1, 2], L[0, 2, 1] = -1, 1 #L_x
    L[1, 0, 2], L[1, 2, 0] = 1, -1 #L_y
    L[2, 0, 1], L[2, 1, 0] = -1, 1 #L_z
    
    Lxhx = L[0, :, :] * h[0]
    Lyhy = L[1, :, :] * h[1]
    Lzhz = L[2, :, :] * h[2]
    return  -(Lxhx + Lyhy + Lzhz)


def reemplazo_cramer(A, b, i):
    '''
    Reemplaza la i-esima columna de A,
    por el vector b.
    '''
    B = np.zeros_like(A)
    B[:,:] = A[:,:]
    B[:, i] = b
    return B


def det_dim2(A):
    return A[0,0]*A[1,1]-A[0,1]*A[1,0]


def det_dim3(A):
    A1 = np.zeros((2, 2))
    A2 = np.zeros((2, 2))
    A3 = np.zeros((2, 2))
    A1[:, 0], A1[:, 1] = A[1:3, 1], A[1:3, 2]
    A2[:, 0], A2[:, 1] = A[1:3, 2], A[1:3, 0]
    A3[:, 0], A3[:, 1] = A[1:3, 0], A[1:3, 1] 
    det1 = det_dim2(A1)
    det2 = det_dim2(A2)
    det3 = det_dim2(A3)
    det = A[0,0]*det1 + A[0,1]*det2+A[0,2]*det3
    return det
    

def esf_to_cart(r, theta, phi):
    s_x =  np.sin(theta) * np.cos(phi)
    s_y = np.sin(theta) * np.sin(phi)
    s_z = np.cos(theta)
    s = [s_x, s_y, s_z] * r
    return np.array(s)


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


def seleccionar_color(sx, sy, sz):
    cos_th = sz / np.sqrt(sx ** 2 + sy ** 2 + sz ** 2)
    if cos_th > 0.7:
        return (0.95, 0, 0)
    elif cos_th > 0.4:
        return (0.8, 0.8, 0.0)
    elif cos_th > 0.125:
        return (0.6, 1, 0.1)
    elif cos_th > -0.125:
        return (0.1, 0.8, 0.4)
    elif cos_th > -0.4:
        return (0.2, 0.7, 0.9)
    elif cos_th > -0.7:
        return (0, 0.00, 1)
    else:
        return (0.5, 0.00, 1)

def plotear_instante_grilla(spines, indice_t, posiciones, fig):
    '''
    Recibe una matriz con spines en formato
    [Nodo_i, intante_tiempo, componente_spin]
    y grafica la grilla en un instante indice_t, 
    (debe ser un entero) graficando los spines como
    flechas de norma 1, que parten en
    las posiciones de los nodos.
    
    Las posiciones se ingresan en un arreglo 2d de formato:
    [nodo_i, componente_pos] y las componentes son
    bidimensionales (nodos viven en un plano).
    
    Se debe ingresar tambien un objeto figure en el que se
    ploteara la grilla de nodos.
    '''
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
        color = seleccionar_color(sx, sy, sz)
        arrow = Arrow3D([x0, xf], [y0, yf], [z0, zf], mutation_scale=10, lw=3, arrowstyle="-|>", color=color)
        ax.add_artist(arrow)
        ax.plot([x0], [y0], [z0], '.', color=(0, 0, 0), markersize=5)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$')



def guardar_spines(spines, posiciones):
    spines_x = spines[:, :, 0]
    spines_y = spines[:, :, 1]
    spines_z = spines[:, :, 2]
    np.savetxt('arreglos/spines_x.dat', spines_x)
    np.savetxt('arreglos/spines_y.dat', spines_y)
    np.savetxt('arreglos/spines_z.dat', spines_z)
    pos_x = posiciones[:, 0]
    pos_y = posiciones[:, 1]
    np.savetxt('arreglos/pos_x.dat', pos_x)
    np.savetxt('arreglos/pos_y.dat', pos_y)

    
def cargar_spines():
    spines_x = np.loadtxt('arreglos/spines_x.dat')
    spines_y = np.loadtxt('arreglos/spines_y.dat')
    spines_z = np.loadtxt('arreglos/spines_x.dat')
    n_nodos, n_t = np.shape(spines_x)
    spines = np.zeros((n_nodos, n_t, 3))
    spines[:, :, 0] = spines_x
    spines[:, :, 1] = spines_y
    spines[:, :, 2] = spines_z

    pos_x = np.loadtxt('arreglos/pos_x.dat')
    pos_y = np.loadtxt('arreglos/pos_y.dat')
    posiciones = np.zeros((n_nodos, 2))
    posiciones[:, 0] = pos_x
    posiciones[:, 1] = pos_y
    return spines, posiciones