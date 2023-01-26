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
K = 1600000
omega_idot = 0
A = 40

#Calcul de variables externes
Wn = np.sqrt(k*d**2/I)
Wd_op = Wn*np.sqrt(1-C_op**2)
Wd_ed = Wn*np.sqrt(1-C_ed**2)
Ceq_op = 2*C_op*Wn
Ceq_ed = 2*C_ed*Wn
sol1_0 = [0,omega_idot]
sol2_0 = [0,omega_idot]


w = 50
ti = 0
tf = 0.5
npts = 1000
t=np.linspace(ti,tf,npts) 
M = A*np.sin(w*t)

def model(y,t,Meq,C,I,K):
        theta, omega = y
        dydt = [omega, (Meq/I)-(C)*omega-(K/I)*theta]
        return dydt

sol_op1 = np.empty(len(t),)
sol_op1[0]=0
sol_op2 = np.empty(len(t),)
sol_op2[0]=0

for i in range(0,len(t)-1):
    dt = [t[i],t[i+1]]
    sol1 = odeint(model,sol1_0,dt,args=(M[i],Ceq_op,I,K))
    sol2 = odeint(model,sol2_0,dt,args=(M[i],Ceq_ed,I,K))
    
    sol1_0 = sol1[1]
    sol2_0 = sol2[1]
    sol_op1[i+1]=sol1[1,0]
    sol_op2[i+1]=sol2[1,0]

sol_op = plt.plot(t, sol_op1, label='PCU Opérationnelle')
sol_ed = plt.plot(t, sol_op2, label='PCU en défaillance')
plt.title("Oscillation de l'aileron")
plt.xlabel('Temps (s)')
plt.ylabel("Angle de l'aileron (rad)")
plt.legend(loc="upper right")
plt.show(sol_ed, sol_op)












