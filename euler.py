#/usr/bin/ipython3 --pylab

#if __name__ == "__main__":

import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
import pylab
import time as tim

#
m1 = 1
m2 = 1
m3 = 1e-16
g = 1

x1=np.array([-0.5,0,0])
x2=np.array([0.5,0,0])
x3=np.array([0,np.sqrt(0.75),0])

v1=np.array([0,-np.sqrt(0.5),0])
v2=np.array([0,np.sqrt(0.5),0])
v3=np.array([-np.sqrt(0.5)*np.sqrt(0.75)*2,0,0])

dt = 2.5e-2
tmax = 100
tout=0.05
time = 0.0
tlast = 0.0

plt.figure(figsize=(8,8))
plt.gca().set_aspect('equal')

plt.scatter (x1[0],x1[1],color='k')
plt.scatter (x2[0],x2[1],color='r')
plt.scatter (x3[0],x3[1],color='b')

plt.xlim ( (-1.5,1.5) )
plt.ylim ( (-1.5,1.5) )
plt.xlabel ('x')
plt.ylabel ('y')
plt.draw()
plt.show(block=False)
    
while (time < tmax):

    v1old=v1
    v2old=v2
    v3old=v3

    v1 = v1 + dt * g * ( \
        + m2 * (x2-x1) / np.linalg.norm(x2-x1)**3 \
        + m3 * (x3-x1) / np.linalg.norm(x3-x1)**3 \
                    )
    v2 = v2 + dt * g * ( \
        + m1 * (x1-x2) / np.linalg.norm(x1-x2)**3 \
        + m3 * (x3-x2) / np.linalg.norm(x3-x2)**3 \
                    )
    v3 = v3 + dt * g * ( \
        + m1 * (x1-x3) / np.linalg.norm(x1-x3)**3 \
        + m2 * (x2-x3) / np.linalg.norm(x2-x3)**3
                    )

    x1 = x1 + v1old * dt
    x2 = x2 + v2old * dt
    x3 = x3 + v3old * dt

    time = time + dt
    
    if (time > tlast):
        tlast = time + tout
        plt.scatter (x1[0],x1[1],color='k',s=1)
        plt.scatter (x2[0],x2[1],color='r',s=1)
        plt.scatter (x3[0],x3[1],color='b',s=1)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)
        

plt.show(block=True)
        
        
        
