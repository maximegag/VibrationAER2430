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


plot_Xmax_od=plt.plot(r, Xmax_op,"b", label='PCU Opérationnelle')
plot_Xmax_ed=plt.plot(r, Xmax_ed,"r", label='PCU En Défaillance')

w50_op=plt.plot(0.197642,1.04043,"b",marker="o")
w100_op=plt.plot(0.395285,1.18389,"b", marker="o")
wn1_op=plt.plot(1,10,"b",marker="o")

w50_ed=plt.plot(0.197642,1.03519,"r",marker="x")
w100_ed=plt.plot(0.395285,1.15395,"r", marker="x")
wn1_ed=plt.plot(1,2,"r",marker="x")

plt.title("Amplitude relative en fonction de r pour les deux modes de fonctionnement")
plt.xlabel("r=ω/ωn")
plt.ylabel("|kX/F0|")
plt.legend(loc="upper right")
plt.show(plot_Xmax_od,plot_Xmax_ed,wn1_op,w100_op,w50_op)




