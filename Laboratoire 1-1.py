#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:51:03 2023

@author: maximegagnon
"""
import numpy as np
import matplotlib.pyplot as plt

#Variable 
k = 40*10**6
I = 25 
d = 0.2
C_op = 0.05
C_ed = 0.25
I = 25
M = 100

#Calcul de variables externes
Wn = np.sqrt(k*d**2/I)
Wd_op = Wn*np.sqrt(1-C_op**2)
Wd_ed = Wn*np.sqrt(1-C_ed**2)
omega_idot = M / I
print(f"Vitesse initiale = {omega_idot}")


# x as delta en fonction du temps
ti = 0
tf = 0.5
npts = 10000
t=np.linspace(ti,tf,npts)   

y_op = (omega_idot/Wd_op*np.sin(Wd_op*t))*np.exp(-C_op*Wn*t)
sol_op=plt.plot(t,y_op,label="PCU Opérationnelle")
y_ed = (omega_idot/Wd_ed*np.sin(Wd_ed*t))*np.exp(-C_ed*Wn*t)
sol_ed=plt.plot(t,y_ed,label="PCU en défaillance")

#Création du tableau
plt.title("Oscillation de l'aileron")
plt.xlabel('Temps (s)')
plt.ylabel("Angle de l'aileron (rad)")
plt.legend(loc="upper right")
plt.show(sol_op,sol_ed)