import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as lin

class Nodo_spin:
    alpha = 0.05
    
    def __init__(self, s, H,  vecinos):
        self.campo_externo = H
        self.vecinos = vecinos
        self.spin = s
        

def h_intercambio(lista_nodos, i):
    nodo = lista_nodos[i]
    vecinos = nodo.vecinos
    n_vec = np.size(vecinos, 0)
    h = np.zeros((3))
    for j in range(n_vec):
        indice_vec = int(vecinos[j, 0])
        s_vec = lista_nodos[indice_vec].spin
        j_vec = vecinos[j, 1]
        h = h + s_vec*j_vec
    return h


def resolver_LL(dt, H_eff, s_0):
    '''
    Recibe un salto de tiempo, el campo magnetico
    efectivo, y un vector de spin inicial s0.
    Retorna la siguiente iteracion del arreglo spin,
    precesandolo en torno al campo.
    '''
    I = np.zeros((3,3)) #cmatriz identidad
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
    H_eff = - H + alpha * cruz(spin, H)
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

def simular_spines(t_array, lista_nodos0):
    '''
    Recibe un arreglo de tiempos, y una lista
    de nodos del instante inicial.
    La funcion retorna una matriz de 2 dimensiones
    conteniendo las grillas de nodos
    en cada uno de los instantes de t_array.
    El formato de retorno es:
    [instante_tiempo, nodo_i]
    '''
    n = np.size(t_array,0)
    lista_nodos = [lista_nodos0]
    for i in range(1, n):
        dt = t_array[i]- t_array[i-1]
        nueva_lista = avanzar_grilla(dt, lista_nodos[i-1])
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
    return  Lxhx + Lyhy + Lzhz


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