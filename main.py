import numpy as np
import os


print("Please type file to read: ")
ftr=input("> ")

print(f"File selected: {ftr}")
print("To continue please press ENTER or CTRL-C to exit.")
input('> ')

initialdata = np.genfromtxt(ftr, encoding='ASCII', dtype=None, names=('Nombre', 'p.m', 
'Tc [K]', 'Pc [bar]', 'omega', 'kappa1','alfa0', 'alfa1', 'beta1', 'beta2', 'Ctv [cm3/gr]'),
skip_header=1, max_rows=3)


