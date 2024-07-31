#------------------------------------------------------------------------------------------------------------
#MODULOS

import pandas as pd
import numpy as np
from scipy.optimize import fsolve

#--------------------------------------------------------------------------------------------------------------
#TRATAMIENTO DE DATOS
# Cargar los datos de la hoja desde el archivo excel
file_path = r'C:\Users\cesar\Desktop\INVESTIGACION\INV-PERSONAL\GAS NATURAL\AGA8-92DC\coef_aga.xlsx'
sheets = pd.read_excel(file_path, sheet_name=['TABLA B1','TABLA B2','TABLA B3'])
tabla_b1 = sheets['TABLA B1']
tabla_b2 = sheets['TABLA B2']
tabla_b3 = sheets['TABLA B3']

# Lista de compuestos a comparar
N_C = ["Methane", "Ethane", "Propane", "iso-Butane", "n-Butane", "iso-Pentane", "n-Pentane", "n-Hexane","Nitrogen","Carbon Dioxide"]

#EXTRACCIÓN DE TABLA 1
a = tabla_b1['a_n'].values
b = tabla_b1['b_n'].values
c = tabla_b1['c_n'].values
k = tabla_b1['k_n'].values
u = tabla_b1['u_n'].values
g = tabla_b1['g_n'].values
q = tabla_b1['q_n'].values
f = tabla_b1['f_n'].values
s = tabla_b1['s_n'].values
w = tabla_b1['w_n'].values

#TRATAMIENTO DE DATOS TABLA 2
# Filtrar las filas donde el 'Compound' está en N_C
df_filtrado = tabla_b2[tabla_b2['Compound'].isin(N_C)]
# Crear un diccionario para mapear los compuestos a índices
compound_to_index2 = {compound: idx for idx, compound in enumerate(N_C)}
# Reindexar el DataFrame según el orden de N_C
df_filtrado = df_filtrado.set_index('Compound').reindex(N_C).reset_index()
# Extraer los datos a partir de la columna 'C'
# Encontrar la posición de la columna 'Compound' y seleccionar las columnas posteriores
columna_compound_index = tabla_b2.columns.get_loc('Compound')
datos_matriz = df_filtrado.iloc[:, columna_compound_index + 1:].values
# Extraer columnas específicas
E_i = datos_matriz[:, 1]  # Primer columna
K_i = datos_matriz[:, 2]  # Segunda columna
G_i = datos_matriz[:, 3]  # Tercera columna
Q_i = datos_matriz[:, 4]  # Cuarta columna
F_i = datos_matriz[:, 5]  # Quinta columna
S_i = datos_matriz[:, 6]  # Sexta columna
W_i = datos_matriz[:, 7]  # Séptima columna


# TRATAMIENTO DE DATOS TABLA 3
# Crear una lista para almacenar los resultados
combinations = []
# Iterar sobre cada nombre en N_C
for compound_i in N_C:
    # Filtrar 'TABLA B3' por la columna 'i' con el valor de compound_i
    filtered_b3_i = tabla_b3[tabla_b3['i'] == compound_i]
    
    # Iterar sobre las filas filtradas de 'TABLA B3'
    for _, row_j in filtered_b3_i.iterrows():
        compound_j = row_j['j']
        
        # Verificar si compound_j está en N_C
        if compound_j in N_C:
            # Obtener los valores posteriores a la columna 'j'
            columna_j = tabla_b3.columns.get_loc('j')
            data_j = row_j.iloc[columna_j + 1:].values
            # Añadir los resultados a la lista
            combinations.append((compound_i, compound_j, *data_j))
# Convertir la lista de combinaciones a un DataFrame para mejor visualización
df_combinations = pd.DataFrame(combinations, columns=['Compound_i', 'Compound_j', 'E*_ij', 'U_ij', 'K_ij', 'G*_ij'])
# Crear matrices para los parámetros binarios
N = len(N_C)
# Inicializar matrices
E_star = np.ones((N, N))
U_ij = np.ones((N, N))
K_ij = np.ones((N, N))
G_star = np.ones((N, N))
# Crear un diccionario para mapear los nombres a índices
compound_to_index = {compound: idx for idx, compound in enumerate(N_C)}
# Asignar valores específicos de combinations
for _, row in df_combinations.iterrows():
    i = compound_to_index[row['Compound_i']]
    j = compound_to_index[row['Compound_j']]
    E_star[i, j] = row['E*_ij']
    U_ij[i, j] = row['U_ij']
    K_ij[i, j] = row['K_ij']
    G_star[i, j] = row['G*_ij']

#------------------------------------------------------------------------------
#FUNCIONES
class AGA8:
    def __init__(self, x):
        self.x = x
        self.N = len(x)
    
    def calcular_F(self, F_i, x):
        F = np.sum(x**2 * F_i)
        return F

    def calcular_Q(self, Q_i, x):
        Q = np.sum(x * Q_i)
        return Q

    def calcular_G(self, G_star, G_i, x):
        G = np.sum(x * G_i)
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                G += 2 * x[i] * x[j] * (G_star[i, j] - 1) * (G_i[i] + G_i[j])
        return G

    def calcular_U5(self, E_i, x, U_ij):
        U5 = np.sum(x * E_i ** (5/2)) ** 2
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                U5 += x[i] * x[j] * (U_ij[i, j] ** 5 - 1) * (E_i[i] * E_i[j]) ** (5/2)
        return U5

    def calcular_C_star(self, n, a, g, q, f, u, G, Q, F, U, T):
        term1 = (G + 1 - g[n]) ** g[n]
        term2 = (Q**2 + 1 - q[n]) ** q[n]
        term3 = (F + 1 - f[n]) ** f[n]
        term4 = U ** u[n]
        term5 = T ** (-u[n])
        C_star_n = a[n] * term1 * term2 * term3 * term4 * term5
        return C_star_n

    def calcular_parametros_binarios(self, E_star, G_star, E_i, G_i):
        E_ij = np.zeros((self.N, self.N))
        G_ij = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                E_ij[i, j] = E_star[i, j] * (E_i[i] * E_i[j]) ** (1/2)
                G_ij[i, j] = G_star[i, j] * (G_i[i] + G_i[j]) / 2
        return E_ij, G_ij

    def calcular_B_star_nij(self, G_ij, Q_i, F_i, S_i, W_i, g, q, f, s, w, n, i, j):
        B_star_nij = (G_ij[i, j] + 1 - g[n]) ** g[n] * (Q_i[i] * Q_i[j] + 1 - q[n]) ** q[n] * \
                    (F_i[i] ** 0.5 * F_i[j] ** 0.5 + 1 - f[n]) ** f[n] * \
                    (S_i[i] * S_i[j] + 1 - s[n]) ** s[n] * (W_i[i] * W_i[j] + 1 - w[n]) ** w[n]
        return B_star_nij

    def calcular_B(self, T, E_ij, K_ij, G_ij, Q_i, F_i, S_i, W_i, a, u, g, q, f, s, w):
        B=0
        for n in range(18):
            sum_ij = 0
            for i in range(self.N):
                for j in range(self.N):
                    B_star_nij = self.calcular_B_star_nij(G_ij, Q_i, F_i, S_i, W_i, g, q, f, s, w, n, i, j)
                    sum_ij += B_star_nij * E_ij[i, j] ** (u[n])  * (K_i[i] * K_i[j]) ** 1.5
            B += a[n] * T ** (-u[n]) * sum_ij
        return B

    def calcular_rho_m(self, P, T, x, R, E_star, G_star, E_i, G_i, a, g, q, f, u, K_i, K_ij, S_i, W_i, s, w):
        def ecuacion_rho_m(rho_m):
            F = self.calcular_F(F_i, x)
            Q = self.calcular_Q(Q_i, x)
            G = self.calcular_G(G_star, G_i, x)
            U = self.calcular_U5(E_i, x, K_ij) ** (1/5)
            E_ij, G_ij = self.calcular_parametros_binarios(E_star, G_star, E_i, G_i)
            K = self.calcular_U5(K_i, x, K_ij) ** (1/5)
            rho_r = K ** 3 * rho_m

            B = self.calcular_B(T, E_ij, K_ij, G_ij, Q_i, F_i, S_i, W_i, a, u, g, q, f, s, w)

            C_star = np.zeros(58)
            for n in range(13, 58):
                C_star[n] = self.calcular_C_star(n, a, g, q, f, u, G, Q, F, U, T)
            
            term_sum = 0
            for n in range(13, 58):
                b_n = b[n]  # Asigna los valores reales de b_n
                c_n = c[n]  # Asigna los valores reales de c_n
                k_n = k[n]  # Asigna los valores reales de k_n
                term_sum += C_star[n] * (b_n - c_n * k_n * rho_r**b_n) * rho_r**b_n * np.exp(-c_n * rho_r**k_n)

            

            Z=(1 + rho_r * B - rho_r * np.sum(C_star[13:19]) + term_sum)
            return R * T * rho_m * Z - P

        rho_m = fsolve(ecuacion_rho_m,0, xtol=1e-8)[0]
        return rho_m

    def calcular_factor_de_compresibilidad_Z(self, P, T, x, R, E_star, G_star, E_i, G_i, a, g, q, f, u, K_i, K_ij, S_i, W_i, s, w):
        rho_m = self.calcular_rho_m(P, T, x, R, E_star, G_star, E_i, G_i, a, g, q, f, u, K_i, K_ij, S_i, W_i, s, w)
        E_ij, G_ij = self.calcular_parametros_binarios(E_star, G_star, E_i, G_i)
        
        U = self.calcular_U5(E_i, x, K_ij) ** (1/5)
        G = self.calcular_G(G_star, G_i, x)
        Q = self.calcular_Q(Q_i, x)
        F = self.calcular_F(F_i, x)
        
        C_star = np.zeros(58)
        for n in range(13, 58):
            C_star[n] = self.calcular_C_star(n, a, g, q, f, u, G, Q, F, U, T)
        
        K = self.calcular_U5(K_i, x, K_ij) ** (1/5)
        rho_r = K ** 3 * rho_m
        
        B = self.calcular_B(T, E_ij, K_ij, G_ij, Q_i, F_i, S_i, W_i, a, u, g, q, f, s, w)
        term_sum = 0
        for n in range(13, 58):
            b_n = b[n]  # Asigna los valores reales de b_n
            c_n = c[n]  # Asigna los valores reales de c_n
            k_n = k[n]  # Asigna los valores reales de k_n
            term_sum += C_star[n] * (b_n - c_n * k_n * rho_r**b_n) * rho_r**b_n * np.exp(-c_n * rho_r**k_n)

        Z = 1 + rho_r * B - rho_r * np.sum(C_star[13:19]) + term_sum
        return Z
    
#------------------------------------------------------------------------------------------------
#PRUEBAS CON EJEMPLOS
# Fracciones molares de los compuestos (ejemplo)
P= 6657160
T= 225
R = 8.314  # Constante de los gases
x = np.array([0.90644,0.04553,0.00833,0.001,0.00156,0.0003,0.00045,0.0004,0.03134,0.00466])
# Crear instancia de la clase
modelo = AGA8(x)
# Calcular el factor de compresibilidad Z
ZAGA = modelo.calcular_factor_de_compresibilidad_Z(P, T, x, R, E_star, G_star, E_i, G_i, a, g, q, f, u, K_i, K_ij, S_i, W_i, s, w)
print("Factor de Compresibilidad AGA8: ",ZAGA)