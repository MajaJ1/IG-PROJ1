# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 19:09:03 2022

@author: admin
"""
import math
import numpy as np

X = 3744087.335 
Y = 954562.958 
Z = 5057780.953

a = 6378137
e2 = 0.00669438002290



fi =  0.9216686397488353
lam = 0.2496340124796006


s = 53000
A = (90+(00/60)+(0/3600))*math.pi/180
z = (90+(00/60)+(0/3600))*math.pi/180
print('n_e_u')
def s_A_z2neu(s, A, z):
    #A = np.deg2rad(A)
    #z = np.deg2rad(z)
    n = s*np.sin(z)*np.cos(A)
    e = s*np.sin(z)*np.sin(A)
    u = s*np.cos(z)
    return n, e, u
n, e, u = s_A_z2neu(s, A, z)
print('n=', n, 'e=', e, 'u=', u)
def neu2dXYZ(n, e, u, fi, lam):
    
    R = np.array([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
              [-np.sin(fi)*np.sin(lam), np.cos(lam), np.cos(fi)*np.sin(lam)],
              [np.cos(fi), 0, np.sin(fi)]])
    dx = np.linalg.inv(R.transpose()) @ np.array([n, e, u])
    return dx

dx = neu2dXYZ(n, e, u, fi, lam)

Xb = X + dx[0]
Yb = Y + dx[1]
Zb = Z + dx[2]

print('Xc = ', Xb,  ' Yb = ' ,Yb, 'Z = ' , Zb)