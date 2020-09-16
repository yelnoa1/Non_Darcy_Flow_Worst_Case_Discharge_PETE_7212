# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 21:53:28 2020

@author: yelnoa1
"""

from PVTOil import mu_O, Rs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_file = pd.read_csv('vlp.csv')

api = 31.5
t = 200
sept = 4
sepp= 2200
rho = 56.8
k = 600
h = 100
re = 4000
rw = 0.34517
p_avg = 12000
pbp = 3400
sg = 0.6
diam = 8.5
#j = 50

def Darcy(q,p_avg):
    temp = 0
    pwfs = p_avg #initial guess
    tolerance = 1e-5
    i = 0
    while np.abs(((pwfs - temp) / pwfs)) > tolerance:
        temp = pwfs 
        p_bar = (p_avg + pwfs)/2.0
        rs = Rs(p_bar, api, t, sg, sept, sepp)
        mu_o = mu_O(api, p_bar, t, rs, pbp)
        F = rs*(sg/api)**0.5 + 1.25*t
        fvf = 0.972 + 0.000147*F**1.175
        pwfs = p_avg -( q * 141.2 * mu_o * fvf * (np.log( 0.472 * re / rw )- 0.75) / (k * h) )
        i+=1
    print(f"Number of iterations = {i}")
    return pwfs
    

def NonDarcy(q):
    temp = 0
    pwf = p_avg #initial guess
    tolerance = 1e-5
    i = 0
    while np.abs(((pwf - temp) / pwf)) > tolerance:
        temp = pwf 
        p_bar = (p_avg + pwf)/2
        rs = Rs(p_bar, api, t, sg, sept, sepp) #assumes equal to GOR
        mu_o = mu_O(api, p_bar, t, rs, pbp)
        F = rs*(sg/api)**0.5 + 1.25*t
        fvf = 0.972 + 0.000147*F**1.175
        c = (mu_o*fvf/(0.001127*2*3.14*k*h))*(np.log(0.472*re/rw))
        beta = (6.15*10**10)/(k**1.55)
        D = (9.08*10**-13*beta*fvf**2*rho)/(4*3.14**2*h**2*rw)
        del_p = q*c+(D*(q**2))
        pwf = p_avg - del_p
        i+=1
    print(f"Number of iterations = {i}")
    return pwf



q = np.arange(0,1100000,30000)
pw=[]
pwf_j = []
for i in q:
    pw.append(NonDarcy(i))
    pwf_j.append(Darcy(i,p_avg))

time = data_file['x']
timetest = time[time<60000]
prod = data_file['y']
prodtest = prod[time<700000]
plt.plot(q,pw,'m-', label='Non-Darcy')
plt.plot(q,pwf_j,'y-', label = 'Darcy')
plt.plot(time,prod,'b-', label = 'VLP')
font = {'family' : 'times new roman','weight' : 'normal','size'   : 12}
plt.rc('font', **font)
plt.ylabel('Bottom hole pressure (psia)', fontsize=12, fontweight = 'normal')
plt.xlabel('Flow rate (bbl/day)', fontsize=12,fontweight = 'normal')
plt.ylim(0,15000)
plt.xlim(0,1000000)
plt.grid(which='major', linestyle=':', linewidth='0.5', color='red')
plt.legend(loc='upper right', borderaxespad=0.)
#plt.savefig('pwfvsq.png',bbox_inches='tight', dpi = (600))
