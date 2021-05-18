# -*- coding: utf-8 -*-
"""
Created on Mon May 10 22:10:09 2021

@author: linch
"""
# this code use solve_ivp to solve complex/stiff ode based on Tunisia model
# the output includes one plot with full scale and one subplot with individual concentration vs time

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# create a complex ODE function of concentration based on time variable
def ode_osci(t,x): # the order depends on solve_ivp
    # constant: rate constant for each reactions based on ref-Tunisia model
    R1 = 1.43 * 10**3  # M-3s-1
    R2 = 2*10**10      # M-2s-1
    R3 = 3.1*10**12    # M-2s-1
    R3r = 2.2          # s-1
    R4 = 7.3*10**3     # M-2s-1
    R4r = 1.7*10**7    # 
    R5 = 6*10**5       # M-1s-1
    R6 = 1*10**4       # M-1s-1
    R7  = 3.2*10**4    # M-1s-1 
    R8  = 7.5*10**5    # M-1s-1
    R9 =  40           # M-1s-1
    R10 = 37           # M-1s-1
    
    # assign each ODE to a vector element
    # initial concentration of bottle A from x array
    A = x[0]  # I-
    B = x[1]  # I2
    C = x[2]  # IO3-            
    D = x[3]  # HIO2         
    E = x[4]  # HOI       
    F = x[5]  # IO2       
    G = x[6]  # MnOH2+                       
    I = x[7]  # HO2
    J = x[8]  # MA 
    K = x[9]   # H2O2     

    # define each ODE
    dA = - R1*(H**2)*A*C - R2*H*A*D - R3*E*H*A + R3r*B + R9*B*J/(1+C9*B) + R10*E*K +k0*(y1o-A)
    dB = R3*E*H*A - R3r*B - R9*B*J/(1+C9*B) + k0*(y2o-B)
    dC = - R1*(H**2)*A*C - R4*C*D*H + R4r*(F*F) + R5*(D*D) + k0*(y3o-C)
    dD = R1*(H**2)*A*C - R2*H*A*D - R4*C*D*H + R4r*(F*F) - 2*R5*(D*D) + R6*F*(MnC - G) + k0*(y4o-D)
    dE = R1*(H**2)*A*C + 2*R2*H*A*D - R3*E*H*A + R3r*B + R5*(D*D) - R10*E*K + k0*(y5o-E)
    dF =  2*R4*C*D*H - 2*R4r*(F*F) - R6*F*(MnC - G) + k0*(y6o-F)
    dG = R6*F*(MnC - G) - R7*G*K + k0*(y7o-G)
    dI = R7*G*K - 2*R8*(I*I) + k0*(y8o-I)
    dJ = -R9*B*J/(1+C9*B) + k0*(y9o-J)
    dK = -R7*G*K + R8*(I*I) - R10*E*K + k0*(y10o-K)
    dydt = np.array([dA, dB, dC, dD, dE, dF, dG, dI, dJ, dK]) # store all the value in an array and return
    return dydt

# initial conditions in bottle A
#default x0 = np.array([10**(-8), 6*10**(-7), 10**(-2), 10**(-10), 10**(-10), 10**(-10), 10**(-13), 0, 10**(-3), 0])
x1o = 10**(-8)      # I-
x2o = 6*10**(-7)    # I2 
x3o = 10**(-2)      # IO3-  
x4o = 10**(-10)     # HIO2 
x5o = 10**(-10)     # HOI
x6o = 10**(-10)     # IO2 
x7o = 10**(-13)     # MnOH2+
x8o = 0             # HO2
x9o = 10**(-3)      # MA   
x10o = 0            # H2O2 
x0 = np.array([x1o, x2o, x3o, x4o, x5o, x6o, x7o, x8o, x9o, x10o]) # store initial condition of bottle A in an array
str_name = ["I-","I2", "IO3-", "HIO2", "HOI", "IO2 ", "MnOH2+", "HO2", "MA", "H2O2"] # store name of each chemicals for future plot

#initial concentration of bottle B which will be added into bottle A, 
#this value will be affected by flow rate
y1o = 0           # I-
y2o = 10**(-6)    # I2      10**(-6)
y3o = 0.035       # IO3-    0.035
y4o = 0           # HIO2 
y5o = 0           # HOI
y6o = 0           # IO2 
y7o = 0           # MnOH2+
y8o = 0           # HO2
y9o = 0.0015      # MA      0.0015
y10o = 0.33       # H2O2    0.33
H = 0.056  # constant concentration of H+ due to no pH changes
k0 = 1/156 # contant flow rate 
C9 = 10**4 # M-1 contant derived from reaction 9
MnC = 0.004  # contant = [Mn2+] + [MnOH2+]

# solve the ODE function and store it in arrays

# declare a time vector (time window)
t_span = np.array([0,5000])  # set time from 0 ~ 5000 sec
times = np.linspace(t_span[0],t_span[1],5001) # set the time_step = 1
# Use Radau method in solve_ivp to solve Stiff systems of ODES
x = solve_ivp(ode_osci, t_span, x0, method = 'Radau', t_eval = times,  rtol = 1e-8 ,atol = 1e-6)
print(x) # print out the x to see if it fail/seccess # can comment or not

# store the solved ODE concentration of each chemicals to the array for further plot
t = x.t
A = x.y[0]
B = x.y[1]
C = x.y[2]
D = x.y[3]
E = x.y[4]
F = x.y[5]
G = x.y[6]
I = x.y[7]
J = x.y[8]
K = x.y[9]

# plot the result
# plot the change of concentration, x axis: time, y axis: concentration; 

# first plot whole concentration in one plot in full scale
plt.figure(figsize = (14,8)) # set the fig size
plt.rcParams['font.size'] = '24' # set the same font size
plt.plot(t,A, label = 'I-') # plot the concentration vs time with assigned label
plt.plot(t,B, label = 'I2')
plt.plot(t,C, label = 'IO3')
plt.plot(t,D, label = 'HIO2')
plt.plot(t,E, label = 'HOI')
plt.plot(t,F, label = 'IO2')
plt.plot(t,G, label = 'MnOH2+')
plt.plot(t,I, label = 'HO2')
plt.plot(t,J, label = 'MA')
plt.plot(t,K, label = 'H2O2')
plt.title("Concentration versus Time") # title
plt.xlabel("time (s)") # x label
plt.ylabel("Concentration (M)") # y label
plt.legend(loc = "upper left") # set the location of label
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # set the label as scientific notation
plt.show()


# second plot: use subplot/loop to plot each individual concentration
fig = plt.figure(figsize = (25,10))  # set larger fig size to contain 10 subplot
plt.rcParams['font.size'] = '12'     # set the same font size 12
for i in range(0,len(str_name)):     # loop for each chemicals
    plt.subplot(2,5,i+1)  # set subplot differnece from "i"
    plt.plot(t,x.y[i], label = str_name[i])  # plot subplot with label
    plt.xlabel("time (s)")
    plt.ylabel("Concentration (M)")
    plt.legend(loc = "upper left") # set the location of label
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # set the label as scientific notation
plt.subplots_adjust(wspace = 0.5, hspace = 0.5) # set the distance (wspace/hspace) between each subplot
plt.suptitle("Concentration versus Time")
plt.show()



