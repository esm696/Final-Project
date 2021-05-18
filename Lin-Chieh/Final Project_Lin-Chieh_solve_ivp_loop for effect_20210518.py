# -*- coding: utf-8 -*-
"""
Created on Mon May 10 22:10:09 2021

@author: linch
"""
# this code use solve_ivp to solve complex/stiff ode based on Tunisia model
# the output includes plots with concentration effect with different factor

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# create a complex ODE function of concentration based on time variable
def ode_osci(t,x): # the order depends on solve_ivp
    # constant: rate constant for each reactions based on ref
    R1 = 1.43 * 10**3  # M-3s-1
    R2 = 2*10**10      # M-2s-1
    R3 = 3.1*10**12    # M-2s-1
    R3r = 2.2          # s-1
    R4 = 7.3*10**3     # M-2s-1
    R4r = 1.7*10**7    # not influence significantly
    R5 = 6*10**5       # M-1s-1
    R6 = 1*10**4       # M-1s-1
    #R6r = 10**4        # the reverse reaction has been neglected
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
x0 = np.array([x1o, x2o, x3o, x4o, x5o, x6o, x7o, x8o, x9o, x10o])

#initial concentration of bottle B as well as H+ concentration and flow rate
y0 = np.array([0, 10**(-6),0.035,0,0,0,0,0,0.0015,0.33, 0.056, 1/156]) 
str_name = ["I-","I2", "IO3-", "HIO2", "HOI", "IO2 ", "MnOH2+", "HO2", "MA", "H2O2", "H+", "flow rate"] # store name of each chemicals for future plot
str_count = len(str_name)  # length of str_name for future loop
H = 0.056  # constant concentration of H+ due to no pH changes 
k0 = 1/156 # contant flow rate 
C9 = 10**4 # M-1 contant derived from reaction 9
MnC = 0.004  # contant = [Mn2+] + [MnOH2+]

# solve the ODE function and loop it for find out the concentration effect
fac = [0.25,0.5,0.75,1,1.25] # the parameter will time fac to find out the effect
count = len(fac)

# this array is the temparary variable which store the y0 and times the fac
# the purpose here is not to change the initial value (constant), and can get the constant initial value for each loop
a = np.zeros(str_count)   

# loop for each parameter y0
for i in range(0,str_count): # loop for each parameter: chemicals concentration and flow rate
    if y0[i] ==0: # if the value, not run the loop since it is unmeaningful
        continue
    
    fig = plt.figure(figsize = (14,8)) # create a fig to contain one plot
    plt.rcParams['font.size'] = '12'   # set the constant font size as 12
    for j in range(0,count): # loop for each fac,
        for k in range(0,str_count):  # loop for assign y0 to a array
            a[k] = y0[k]
            if k == i:   # if k = i, then this value will time fac
                a[k] = y0[k]*fac[j]
        # after the loop, there is a new a array which will time fac value for one parameter
        # assign the a array back to yio variables, and use those varialbes to solve ODE
        y1o = a[0]    # I-     0
        y2o = a[1]    # I2     %%
        y3o = a[2]    # IO3-   %%
        y4o = a[3]    # HIO2   0
        y5o = a[4]    # HOI    0
        y6o = a[5]    # IO2    0
        y7o = a[6]    # MnOH2+ 0
        y8o = a[7]    # HO2    0
        y9o = a[8]    # MA     %%
        y10o = a[9]   # H2O2   %%
        H = a[10]     # H+     %%
        k0 = a[11]    # flow rate%%
        
        # declare a time vector (time window)
        t_span = np.array([0,5000])  # set time from 0 ~ 5000 sec
        times = np.linspace(t_span[0],t_span[1],5001) # set the time_step = 1
        # Use Radau method in solve_ivp to solve Stiff systems of ODES
        x = solve_ivp(ode_osci, t_span, x0, method = 'Radau', t_eval = times,  rtol = 1e-8 ,atol = 1e-6)
        t = x.t # time array
        # plot concentration of I- for each fac during the loop
        plt.plot(t,x.y[0], label = 'I- with ' + str(fac[j])+ '*' + str_name[i])
    plt.xlabel("time (s)")
    plt.ylabel("Concentration (M)")
    plt.legend(loc = "upper left")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # set the label as scientific notation
    plt.title("Concentration versus Time")
    plt.show()