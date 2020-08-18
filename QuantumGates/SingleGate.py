#!/usr/bin/env python
# coding: utf-8

# In[1]:


from qutip import *


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np


# In[3]:


# Dot product between a matrix vector and a vector  
def dot(x,y):
    result = 0

    for i in range(len(y)):
        result = result + x[i]*y[i]

    return result


# In[4]:


#Hamiltonian
def hamiltonian(hbar,gyro_ratio,B,amp):
        
    # Magnetic fields direction
    B = np.array([1,0,1])

    # Magnetic field
    for i in range(len(B)):
        B[i] = B[i]*amp[i]

    # Operators

    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()

    # Spin angular momentum operator
    I = [sx,sy,sz]

    # Magnetic moment operator
    mu = I

    for i  in range(len(I)):
        mu[i] = hbar * gyro_ratio * I[i]

    # Hamiltonian
    H = -1*dot(mu,B)
        
    return H


# In[5]:


class SingleGate:
    
    def __init__(self, hbar, gyro_ratio, B = np.array([1,0,1]), amp=[1,1,1]):
        
        self.H = hamiltonian(hbar,gyro_ratio,B,amp)
        self.gyro_ratio = gyro_ratio
        self.B = B
        self.amp = amp
    
    
    # Pulse and amplitude
    def pulse(self, t, psi):

        fac = 10 
        n = fac * t
        tlist = np.linspace(0, t, n)
        res = mesolve(self.H, psi, tlist, [], [])

        return res, tlist, n
    
    def measure(self, qbit, n_measurements, true=0, false=0):
        a = qbit[0]  # a
        b = qbit[0]  # b
        
        bins = []
        n_false = 0
        error = 0
        
        for i in range(n_measurements):
            
            if np.random.random() <= np.abs(a) ** 2:
                result = 0
            else:
                result = 1

            bit = result
            if bit == false:
                n_false += 1

            bins.append(bit)

        error = n_false/n_measurements
    
        return bins, error


# In[6]:


#Probing the results
def probe(res, t_probe, t, n):
    n_probe = (t_probe/t)*(n)
    if (t_probe/t)==1 :
        return res.states[-1]
    else:
        return res.states[int(n_probe)]


# In[7]:


# Making a direction from qubit
def qvec(vector):
    
    theta = []
    phi = []

    for i in range(len(vector)):
        a, b = vector[i][0],vector[i][1]

        theta_i = 2 * np.arccos(np.absolute(a))

        if theta_i!=0.0 and b!=1.0 : 
            phi_i = np.arccos((b + np.conj(b))/((2.0)*np.sin(theta_i/2)))
            phi_i = phi_i + np.arccos((a + np.conj(a))/((2.0)*np.cos(theta_i/2)))
        else:
            phi_i = 0.0

        theta.append(theta_i)
        phi.append(np.absolute(phi_i))
            
    theta, phi = np.asarray(theta), np.asarray(phi)
    theta = theta.flatten()
    phi = phi.flatten()
        
    return theta, phi


# In[8]:


# Direction of any Vector
def direction(vector): 
    x, y, z = vector[0],vector[1],vector[2]
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)
    theta = np.arccos(z / r)
    return np.array([theta, phi])


# In[9]:


# Making a Qubit
def qubit(angle):
    up = ket("0")
    down = ket("1")
    a = np.cos(angle[0]/2.0)
    b = np.exp(1j*angle[1])*np.sin(angle[0]/2.0)
    qubit = a*up + b*down
    
    return qubit

