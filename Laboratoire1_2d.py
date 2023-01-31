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

#Calcul de la fréquence naturelle
Wn = np.sqrt(k*d**2/I)
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

#Calcul de variables externes
r=np.linspace(ti,tf,npts)  

Xmax_op=1/(np.sqrt((1-(r)**2)**2+(2*C_op*(r))**2))
Xmax_ed=1/(np.sqrt((1-(r)**2)**2+(2*C_ed*(r))**2))


plot_Xmax_od=plt.plot(r, Xmax_op,"b", label='PCU Opérationnelle')
plot_Xmax_ed=plt.plot(r, Xmax_ed,"r", label='PCU En Défaillance')

w50_op=plt.plot(r50,A50_op,"b",marker="o")
w100_op=plt.plot(r100,A100_op,"b", marker="o")
wn1_op=plt.plot(rn1,An1_op,"b",marker="o")

w50_ed=plt.plot(r50,A50_ed,"r",marker="x")
w100_ed=plt.plot(r100,A100_ed,"r", marker="x")
wn1_ed=plt.plot(rn1,An1_ed,"r",marker="x")

#plt.title("Amplitude relative en fonction du ratio de fréquen")
plt.xlabel("r=ω/ωn")
plt.ylabel("|kX/F0|")
plt.legend(loc="upper right")
plt.show()




