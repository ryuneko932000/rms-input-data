import numpy as np
import os


print("Escriba el archivo para leer: ")
ftr=input("> ")

print(f"Archivo seleccionado: {ftr}")
print("Para continuar presione ENTER o CTRL-C para salir.")
input('> ')

data = np.genfromtxt(ftr, encoding='UTF-8',dtype=None, names=('nombre', 'pm', 
'tc', 'pc', 'omega', 'kappa1','alfa0', 'alfa1', 'beta1', 'beta2', 'ctv'),
skip_header=1, max_rows=1)


compuesto = data['nombre']
pm = data['pm']
tc = data['tc']
pc = data['pc']
omega = data['omega']
kappa1 = data['kappa1']
alfa0 = data['alfa0']
alfa1 = data['alfa1']
beta1 = data['beta1']
beta2 = data['beta2']
ctv = data['ctv']

print(compuesto, pm, tc, pc, omega, kappa1, alfa0, alfa1, beta1, beta2, ctv)








