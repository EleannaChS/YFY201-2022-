# SqBilliardQM.py
# Eleanna Choraiti Sideri
# AEM 4406

import numpy as np
import matplotlib.pyplot as plt

Nmax = 101
dt = 0.01
dx = 0.2
dx2 = dx*dx 
fc = dt/dx2
Tmax = 200
Kx = 10.
Ky = 15. 
xin = 50.
yin = 50.

I = np.zeros((Nmax,Nmax),float) 
R = np.zeros((Nmax,Nmax),float)
for i in range(1,Nmax-1): #Initial psi
    for j in range(1 ,Nmax-1):
        Gauss = np.exp(-.05*(i-yin)**2-.05*(j-xin)**2)
        R[i,j] = Gauss*np.cos(Kx*j + Ky*i)
        I[i,j] = Gauss*np.sin(Kx*j + Ky*i)
for t in range (0, Tmax): # Step through time
    R[1:-1 ,1:-1] = R[1:-1 ,1:-1] - fc*(I[2:,1:-1] + I[0:-2 ,1:-1]
                   -4*I[1:-1,1:-1] + I[1:-1,2:] + I[1:-1 ,0:-2])
    I[1:-1 ,1:-1] = I[1:-1 ,1:-1] + fc*(R[2:,1:-1] + R[0:-2,1:-1]
                    -4*R[1:-1,1:-1] + R[1:-1,2:] + R[1:-1,0:-2])
x = y = np.arange(0, Nmax)
X, Y = np.meshgrid (x,y)
Z = (I**2 + R**2)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection="3d")
ax.plot_wireframe(X ,Y ,Z, color='crimson')
ax.set_xlabel("™x")
ax.set_ylabel("™y")
ax.set_zlabel("$\Psi$")
ax.set_title("$\Psi$ (x ,y, t ) at 200 Times on Square Billiard Table")
plt.show()