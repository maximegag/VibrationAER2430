#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:51:03 2023

@author: maximegagnon
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
I = 25
M = 0
K = 1600000

#Calcul de variables externes
Wn = np.sqrt(k*d**2/I)
Wd_op = Wn*np.sqrt(1-C_op**2)
Wd_ed = Wn*np.sqrt(1-C_ed**2)



Ceq_op = 2*C_op*Wn
Ceq_ed = 2*C_ed*Wn
omega_idot = 4

# x as delta en fonction du temps
ti = 0
tf = 0.5
npts = 10000



def model(y,t,M,C,I,K):
    theta, omega = y
    dydt = [omega, (M/I)-(C)*omega-(K/I)*theta]
    
    return dydt

y0 = [0,omega_idot]
t=np.linspace(ti,tf,npts)   

sol_op = odeint(model,y0,t,args=(M,Ceq_op,I,K))
sol_ed= odeint(model,y0,t,args=(M,Ceq_ed,I,K))


sol_op = plt.plot(t, sol_op[:, 0], label='PCU Opérationnelle')
sol_ed = plt.plot(t, sol_ed[:, 0], label='PCU en défaillance')
plt.title("Oscillation de l'aileron")
plt.xlabel('Temps (s)')
plt.ylabel("Angle de l'aileron (rad)")
plt.legend(loc="upper right")
plt.show(sol_ed,sol_op)













