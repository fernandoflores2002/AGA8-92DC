#-----------------------------------------------------------------------------------------
#MODULOS
import math
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------------------
#TRATAMIENTO DE DATOS

import pandas as pd
file_path = 'C:/Users/cesar/Desktop/Proyecto-Gas Natural/coef_aga.xlsx'
sheets = pd.read_excel(file_path, sheet_name=['TABLA B1','TABLA B2','TABLA B3'])
tabla_b1 = sheets['TABLA B1']
tabla_b2 = sheets['TABLA B2']
tabla_b3 = sheets['TABLA B3']

#TABLA B1 (la tabla 1 se puede evaluar puesto que no depende del numero de componentes)
#De 0 a 57 (cada arreglo)
a_n = tabla_b1['a_n'].values.tolist()
b_n = tabla_b1['b_n'].values.tolist()
c_n = tabla_b1['c_n'].values.tolist()
k_n = tabla_b1['k_n'].values.tolist()
u_n = tabla_b1['u_n'].values.tolist()
g_n = tabla_b1['g_n'].values.tolist()
q_n = tabla_b1['q_n'].values.tolist()
f_n = tabla_b1['f_n'].values.tolist()
s_n = tabla_b1['s_n'].values.tolist()
w_n = tabla_b1['w_n'].values.tolist()

#-----------------------------------------------------------------------------------------
#BASE DE DATOS PARA LOS 10 COMPONENTES DE LA COMPOSICION
#x=["Methane","Ethane","Propane","iso-Butane","neo-Butane","iso-Pentane","neo-Pentane","neo-Hexane","Nitrogen","Carbon Dioxide"]
#x=[1,4,5,11,12,13,14,15,2,3]
N_C = ["Methane", "Ethane", "Propane", "iso-Butane", "n-Butane", "iso-Pentane", "n-Pentane", "n-Hexane","Nitrogen","Carbon Dioxide"]
#x=[0.9658,0.01815,0.00405,0.00099,0.00102,0.00047,0.00032,0.00063,0.00269,0.00589]
#T=275 #K
#P=2175940 #Pa
R=8314.5 # Pa*m3/mol*K

#TABLA B2
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

#TABLA B3
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
Estar_ij = np.ones((N, N))
U_ij = np.ones((N, N))
K_ij = np.ones((N, N))
Gstar_ij = np.ones((N, N))
# Crear un diccionario para mapear los nombres a índices
compound_to_index = {compound: idx for idx, compound in enumerate(N_C)}
# Asignar valores específicos de combinations
for _, row in df_combinations.iterrows():
    i = compound_to_index[row['Compound_i']]
    j = compound_to_index[row['Compound_j']]
    Estar_ij[i, j] = row['E*_ij']
    U_ij[i, j] = row['U_ij']
    K_ij[i, j] = row['K_ij']
    Gstar_ij[i, j] = row['G*_ij']

#-----------------------------------------------------------------------------------------

#FUNCIONES
def F(x,F_i):
    F = 0
    for i in range(0,len(x)): 
        F+=x[i]**(2)*F_i[i]
    return F

def Q(x,Q_i):
    Q=0
    for i in range(0,len(x)):
        Q+=x[i]*Q_i[i]
    return Q

def G(x,G_i,Gstar_ij):
    G1=0
    for i in range(0,len(x)):
        G1+=x[i]*G_i[i]
    G2=0
    for i in range(0,len(x)):
        for j in range(i+1,len(x)):
            G2+=x[i]*x[j]*(Gstar_ij[i][j]-1)*(G_i[i]+G_i[j])
    G=G1+2*G2
    return G

def U(x,E_i,U_ij):
    U1=0
    for i in range(0,len(x)):
        U1+=x[i]*E_i[i]**(5/2)
    U1=U1**2
    U2=0
    for i in range(0,len(x)):
        for j in range(i+1,len(x)):
            U2+=x[i]*x[j]*(U_ij[i][j]**5-1)*(E_i[i]+E_i[j])**(5/2)
    U=U1+2*U2
    U=U**(1/5)
    return U

def K(x,K_i,K_ij):
    K1=0
    for i in range(0,len(x)):
        K1+=x[i]*K_i[i]**(5/2)
    K1=K1**2
    K2=0
    for i in range(0,len(x)):
        for j in range(i+1,len(x)):
            K2+=x[i]*x[j]*(K_ij[i][j]**5-1)*(K_i[i]+K_i[j])**(5/2)
    K=K1+2*K2
    K=K**(1/5)
    return K

#LISTA
def Cstar(x,F_i,Q_i,G_i,Gstar_ij,E_i,U_ij,a_n,g_n,q_n,f_n,u_n,T):
    #Variables de inicializacion
    Gf=0
    Qf=0
    Ff=0
    Uf=0
    Cstar=[0]*len(a_n)
    #Calculos
    Gf=G(x,G_i,Gstar_ij)
    Qf=Q(x,Q_i)
    Ff=F(x,F_i)
    Uf=U(x,E_i,U_ij)
    for n in range(0,len(a_n)):
        Cstar[n]=a_n[n]*(Gf+1-g_n[n])**(g_n[n])*(Qf**2+1-q_n[n])**(q_n[n])*(Ff+1-f_n[n])**(f_n[n])*Uf**(u_n[n])*T**((-1)*u_n[n])
    return Cstar

def E_ij(x,Estar_ij,E_i):
    Eij=np.ones((len(x), len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(x)):
            Eij[i][j]=Estar_ij[i][j]*(E_i[i]*E_i[j])**(1/2)
    return Eij


def G_ij(x,Gstar_ij,G_i):
    Gij=np.ones((len(x), len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(x)):
            Gij[i][j]=Gstar_ij[i][j]*(G_i[i]*G_i[j])/2
    return Gij    

#LISTA TRIDIMENSIONAL
def Bstar(x, Gstar_ij, G_i, Q_i, F_i, S_i, W_i, g_n, q_n, f_n, s_n, w_n):
    # Inicializar Bstar con ceros usando listas anidadas
    Bstar = [[[0 for _ in range(len(x))] for _ in range(len(x))] for _ in range(len(g_n))]
    
    # Calcular G_ij
    G_ijf = G_ij(x, Gstar_ij, G_i)
    
    # Llenar Bstar con los valores calculados
    for n in range(0,len(g_n)):
        for i in range(0,len(x)):
            for j in range(0,len(x)):
                Bstar[n][i][j] = (G_ijf[i][j] + 1 - g_n[n]) ** g_n[n] * \
                                 (Q_i[i] * Q_i[j] + 1 - q_n[n]) ** q_n[n] * \
                                 (F_i[i] ** 0.5 * F_i[j] ** 0.5 + 1 - f_n[n]) ** f_n[n] * \
                                 (S_i[i] * S_i[j] + 1 - s_n[n]) ** s_n[n] * \
                                 (W_i[i] * W_i[j] + 1 - w_n[n]) ** w_n[n]
    return Bstar 

def B(x, Gstar_ij, G_i,Estar_ij,E_i, Q_i, F_i, S_i, W_i,a_n, g_n, q_n, f_n, s_n, w_n,u_n,T):
    B=0
    B2=0
    Bstarf=Bstar(x, Gstar_ij, G_i, Q_i, F_i, S_i, W_i, g_n, q_n, f_n, s_n, w_n)
    E_ijf=E_ij(x,Estar_ij,E_i)
    for n in range(0,18):
        B+=a_n[n]*T**(u_n[n])*B2
        for i in range(0,len(x)):
            for j in range(0,len(x)):
                B2+=x[i]*x[j]*Bstarf[n][i][j]*E_ijf[i][j]**(u_n[n])*(K_i[i]*K_i[j])**(3/2)
    return B

def calc_rho_m(x,Gstar_ij,U_ij,K_ij,Estar_ij,G_i,E_i,K_i,Q_i, F_i, S_i, W_i,a_n,b_n,c_n,k_n, g_n, q_n, f_n, s_n, w_n,u_n,T,R,P):
    rho_m=sp.symbols('x')
    Cstarf=Cstar(x,F_i,Q_i,G_i,Gstar_ij,E_i,U_ij,a_n,g_n,q_n,f_n,u_n,T)
    Bf=B(x,Gstar_ij,G_i,Estar_ij,E_i,Q_i,F_i,S_i,W_i,a_n,g_n,q_n,f_n,s_n,w_n,u_n,T)
    Kf=K(x,K_i,K_ij)**3
    interval1=0.55/Kf
    interval2=0.80/Kf
    sum1=0
    for n in range(12,18):
        sum1+=Cstarf[n]
    sum2=0
    for n in range(12,58):
        sum2+=Cstarf[n]*(b_n[n]-c_n[n]*k_n[n]*(Kf*rho_m)**(k_n[n]))*(Kf*rho_m)**(b_n[n])*(sp.exp((-1)*c_n[n]*(Kf*rho_m)**(k_n[n])))
    
    #Funcion final en funcion de la densidad molar
    f_final=rho_m*R*T*(1+Bf*rho_m-(Kf*rho_m)*sum1+sum2)-P
    df_final=sp.diff(f_final)
    # Convertir la función de sympy a numpy
    #f_np = sp.lambdify(rho_m, f_final, 'numpy')
    # Crear un rango de valores para rho_m
    #rho_m_vals = np.linspace(0,20, 400)
    # Evaluar la función en el rango de valores
    #f_vals = f_np(rho_m_vals)
        # Graficar la función con matplotlib
    #plt.plot(rho_m_vals, f_vals, label='$f_{final}(\\rho_m)$')
    #plt.xlabel('$\\rho_m$')
    #plt.ylabel('$f_{final}$')
    #plt.title('Gráfica de $f_{final}(\\rho_m)$')
    #plt.legend()
    #plt.grid(True)
    #plt.show()

    #METODO DE NEWTON-RAPSHON
    # Método de Newton-Raphson
    f = sp.lambdify(rho_m, f_final)
    df = sp.lambdify(rho_m, df_final)
    x0 = interval1
    tol = 0.02
    n = 50  # número máximo de iteraciones
    for k in range(n):
        x1=x0-f(x0)/df(x0)
        if(abs(x1-x0)<tol):
            return x1
        x0=x1

# Leer base de datos
basedatos = pd.read_excel(r'C:\Users\cesar\Desktop\Proyecto-Gas Natural\data_isotermas.xlsx', sheet_name='10comp')
# Lista para almacenar resultados
results = []
# Número de filas
lon = len(basedatos.index)

for i in range(0, lon):
    # Datos del data frame

    T = basedatos.iloc[i, 0]
    P = basedatos.iloc[i, 1]
    x = np.array(basedatos.iloc[i, 2:12].tolist())  # Obtener las composiciones desde la columna 2 hasta 11
    Z_exp= basedatos.iloc[i, 13]

    rho_m=calc_rho_m(x,Gstar_ij,U_ij,K_ij,Estar_ij,G_i,E_i,K_i,Q_i, F_i, S_i, W_i,a_n,b_n,c_n,k_n, g_n, q_n, f_n, s_n, w_n,u_n,T,R,P)
    rho_m=rho_m*1000
    ZAGA=P/(rho_m*R*T)
    # Almacenar los resultados
    results.append([P,T,Z_exp,ZAGA])

# Crear un DataFrame con los resultados
columns = ['P','T','Z(exp)','Z(AGA)']
data = pd.DataFrame(results, columns=columns)

# Mostrar los datos filtrados
print(data)