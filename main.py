import numpy as np
import os


print("Escriba el archivo para leer: ")
ftr=input("> ")

print(f"Archivo seleccionado: {ftr}")
print("Para continuar presione ENTER o CTRL-C para salir.")
input('> ')

initial_data = np.genfromtxt(ftr, encoding='UTF-8',dtype=None, names=('nombre', 'pm', 
'tc', 'pc', 'omega', 'kappa1','alfa0', 'alfa1', 'beta1', 'beta2', 'ctv'),
skip_header=1, max_rows=1)

data2 = np.genfromtxt(ftr, encoding='UTF-8', dtype=None, names=('opcion', 'ede', 'funcionAlfa'), 
skip_header=5, max_rows=7)


compuesto = initial_data['nombre']
pm = initial_data['pm']
tc = initial_data['tc']
pc = initial_data['pc']
omega = initial_data['omega']
kappa1 = initial_data['kappa1']
alfa0 = initial_data['alfa0']
alfa1 = initial_data['alfa1']
beta1 = initial_data['beta1']
beta2 = initial_data['beta2']
ctv = initial_data['ctv']

mtflag_option1 = data2[0]
mtflag_option2 = data2[1]
mtflag_option3 = data2[2]
mtflag_option4 = data2[3]
mtflag_option5 = data2[4]
mtflag_option6 = data2[5]
mtflag_option7 = data2[6]

print(mtflag_option5['ede'])




print(compuesto, pm, tc, pc, omega, kappa1, alfa0, alfa1, beta1, beta2, ctv)










