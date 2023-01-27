# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Import
import scipy.integrate as integrate
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#Variable 
k = 40*10**6
I = 25 
d = 0.2
C_op = 0.05
C_ed = 0.25
A = 40
ti = 0.00000001
tf = 3
npts = 10000

#Calcul de variables externes
Wn = np.sqrt(k*d**2/I)

r=np.linspace(ti,tf,npts)  


Xmax_op=1/(np.sqrt((1-(r)**2)**2+(2*C_op*(r))**2))
Xmax_ed=1/(np.sqrt((1-(r)**2)**2+(2*C_ed*(r))**2))

r_xmax=np.sqrt(1-2*r**2)
rxmax_ed=A/(k*np.sqrt((1-(r_xmax)**2)**2+(2*C_ed*(r_xmax))**2))

plot_Xmax_od=plt.plot(r, Xmax_op, label='PCU Opérationnelle')
plot_Xmax_ed=plt.plot(r, Xmax_ed, label='PCU En Défaillance')
plot_r_xmax=plt.plot(r, rxmax_ed, label='r(Xmax)')
plt.title("Amplitude de l'Oscillation")
plt.xlabel("ratio")
plt.ylabel("Amplitude de l'Oscillation")
plt.legend(loc="upper right")
plt.show(plot_Xmax_od,plot_Xmax_ed,plot_r_xmax)




