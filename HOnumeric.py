# HOnumeric.py
# Eleanna Choraiti Sideri
# AEM 4406

import numpy as np 
import matplotlib.pylab as plt

def rk4Algor(t,h,N,y,f):
    k1=np.zeros(N); k2=np.zeros(N); k3=np.zeros(N); k4=np.zeros(N)
    k1 = h*f(t,y)
    k2 = h*f(t+h/2.,y+k1/2.)
    k3 = h*f(t+h/2.,y+k2/2.)
    k4 = h*f(t+h,y+k3)
    y=y+(k1+2*(k2+k3)+k4)/6.
    return y

rVec = np.zeros((1000),float) # x values for plot
psiVec = np.zeros((1000),float) # Wave function values

fVec = np.zeros(2) #the functions vector
y = np.zeros((2)) #initial conditions vector

n = 5 # n of the bound state energy

def f(x,y): # ODE RHS
    fVec[0] = y[1]
    fVec[1] = -(2*n+1-x**2)*y[0]
    return fVec

if(n%2==0): y[0]=1e-8  # Set parity
else: y[0]=-1e-8
        
y[1] = 1.

i=0
f (0.0, y)  # RHS at r = 0
dr = 0.01

for r in np.arange(-5 ,5 , dr):     # Compute WF steps of dr
    rVec[i] = r
    y = rk4Algor(r ,dr ,2 ,y ,f)
    psiVec[i] = y[0]
    i = i+1       # Advance i 
    
plt.figure()
plt.plot(rVec, psiVec)
plt.grid()
plt.title("Harmonic Oscillator Wave Function n = 5")
plt.xlabel("x")
plt.ylabel("$\psi(x)$")
plt.show()    