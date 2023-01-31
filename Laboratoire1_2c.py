# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 23:45:24 2023

@author: Thomas
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
sol3_0 = [0,omega_idot]
sol4_0 = [0,omega_idot]
sol5_0 = [0,omega_idot]
sol6_0 = [0,omega_idot]


w50 = 50
w100 = 100
wn1 = Wn

#Calculs des amplitudes
r50=w50/Wn
r100=w100/Wn
rn1=wn1/Wn

A50_op=1/(np.sqrt(((1-r50**2)**2)+(2*C_op*r50)**2))
A100_op=1/(np.sqrt(((1-r100**2)**2)+(2*C_op*r100)**2))
An1_op=1/(np.sqrt(((1-rn1**2)**2)+(2*C_op*rn1)**2))

A50_ed=1/(np.sqrt(((1-r50**2)**2)+(2*C_ed*r50)**2))
A100_ed=1/(np.sqrt(((1-r100**2)**2)+(2*C_ed*r100)**2))
An1_ed=1/(np.sqrt(((1-rn1**2)**2)+(2*C_ed*rn1)**2))

ti = 0
tf = 0.5
npts = 1000
t=np.linspace(ti,tf,npts) 
M50 = A*np.sin(w50*t)
M100 = A*np.sin(w100*t)
Mwn = A*np.sin(wn1*t)

def model(y,t,Meq,C,I,K):
        theta, omega = y
        dydt = [omega, (Meq/I)-(C)*omega-(K/I)*theta]
        return dydt

sol_op1 = np.empty(len(t),)
sol_op1[0]=0
sol_op2 = np.empty(len(t),)
sol_op2[0]=0
sol_op3 = np.empty(len(t),)
sol_op3[0]=0

sol_ed1 = np.empty(len(t),)
sol_ed1[0]=0
sol_ed2 = np.empty(len(t),)
sol_ed2[0]=0
sol_ed3 = np.empty(len(t),)
sol_ed3[0]=0

for i in range(0,len(t)-1):
    dt = [t[i],t[i+1]]
    sol1 = odeint(model,sol1_0,dt,args=(M50[i],Ceq_op,I,K))
    sol2 = odeint(model,sol2_0,dt,args=(M100[i],Ceq_op,I,K))
    sol3 = odeint(model,sol3_0,dt,args=(Mwn[i],Ceq_op,I,K))
    sol4 = odeint(model,sol4_0,dt,args=(M50[i],Ceq_ed,I,K))
    sol5 = odeint(model,sol5_0,dt,args=(M100[i],Ceq_ed,I,K))
    sol6 = odeint(model,sol6_0,dt,args=(Mwn[i],Ceq_ed,I,K))
    
    sol1_0 = sol1[1]
    sol2_0 = sol2[1]
    sol3_0 = sol3[1]
    sol_op1[i+1]=sol1[1,0]
    sol_op2[i+1]=sol2[1,0]
    sol_op3[i+1]=sol3[1,0]
    
    sol4_0 = sol4[1]
    sol5_0 = sol5[1]
    sol6_0 = sol6[1]
    sol_ed1[i+1]=sol4[1,0]
    sol_ed2[i+1]=sol5[1,0]
    sol_ed3[i+1]=sol6[1,0]


fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(nrows=6, sharex='col', sharey='row')
ax1.plot(t,sol_op1,"g",label="ζ_op avec w = 50.00 rad/s")
ax3.plot(t,sol_op2,"b",label="ζ_op avec w = 100.00 rad/s")
ax5.plot(t,sol_op3,"r",label="ζ_op avec w = 252.98 rad/s")
ax2.plot(t,sol_ed1,"y",label="ζ_ed avec w = 50.00 rad/s")
ax4.plot(t,sol_ed2,"c",label="ζ_ed avec w = 100.00 rad/s")
ax6.plot(t,sol_ed3,"m",label="ζ_ed avec w = 252.98 rad/s")
plt.xlabel('Temps (s)')
plt.ylabel("                                                                                                       Angle de rotation (rad)\n")
fig.suptitle("Réponse du sytème aux 3 cas d'impulsion pour les 2 modes de fonctionnement", fontsize=12)
fig.legend()
fig.tight_layout()
plt.show