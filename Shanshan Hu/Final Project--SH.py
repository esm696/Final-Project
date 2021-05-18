#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：Final Project--SH
@File ：Final Project--SH.py
@Author ：Shanshan Hu
@Date ：5/18/2021 2:36 PM
'''
import sympy as sym
import numpy as np
######## simulation of TED TMD TSD for AlN
###### set parameters of the material
ax=0.00031111       #### Lattice parameter
cx=0.00049788
v=0.245             #Poisson's coefficient
u=61.2435    #absorption coefficient(cm-1)
####### input parameters
b1=eval(input("Please input the burgers vector for TED: (choose among 0 and ax)"))  #TED burgers vector
b2=eval(input("Please input the burgers vector for TSD: (choose among 0, cx, -cx)")) #TSD burgers vector
h,k,l=1,1,4 #diffraction plane
###### set initial geometry condition
D=3*10**5           #specimen-film distance (um), 30 cm
wave=0.0001267737      #wavelength(um)
d0=1/np.sqrt(4/3*(h**2+h*k+k**2)/(ax**2)+l**2/(cx**2)) # d spacing before deformation
phi=np.arctan((cx/l)/(ax/2))  #angle between diffraction plane and 0001 plane, unit:rad
theta0= 0.1*np.pi/180              #off-axis angle,unit: rad
thetab=np.arcsin(wave/(2*d0))    #bragg angle in perfect crystal,unit:rad
beta=thetab-phi               #angle between incident beam and y[11-20],unit:rad
beta1=np.pi/2-thetab-phi-theta0      #angle between outgoing beam and z axis,unit:rad
p=1/u/(1/np.sin(beta)+1/np.sin(np.pi/2-beta1))*10000     #penetration depth calculation
q0=np.array([-np.cos(beta-theta0),0,-np.sin(beta-theta0)])   # incident unit vector
n=np.array([np.sin(phi+theta0),0,np.cos(phi+theta0)])     # plane normal before distortion
####### set size and resolusion of the simulated result
l1,l2,l3,l4=-25,25,-25,25         #sample size (lengh):l (um)
w=502   #determines the resolution of the image
w2=51   #penetration depth-1 or steps in z direction
x=np.linspace(l1,l2,w)
y=np.linspace(l1,l2,w)                   # size and resolution
z=np.linspace(-p,0,w2)
####### calculation process
X,Y,Z=x.copy(),y.copy(),z.copy()
X,Y,Z=np.meshgrid(X,Y,Z)
x,y,z=sym.symbols("x,y,z")
theta=-theta0                  #inclined angle
bs=b2
x_prime2,y_prime2,z_prime2=-y,x,z

y3=y_prime2*np.cos(theta)-z_prime2*np.sin(theta)
z3=y_prime2*np.sin(theta)+z_prime2*np.cos(theta)      #use this in surface relaxation
y4=-y_prime2*np.cos(theta)-z_prime2*np.sin(theta)
z4=y_prime2*np.sin(theta)-z_prime2*np.cos(theta)

R=sym.sqrt(np.square(x_prime2)+np.square(y_prime2)+np.square(z_prime2))
omega=sym.atan(y_prime2/x_prime2)-sym.atan(y3/x_prime2)+sym.atan(x_prime2*R*np.sin(theta)/(y_prime2*y3+np.square(x_prime2)*np.cos(theta)))
omega1=sym.atan(y_prime2/x_prime2)-sym.atan(y4/x_prime2)+sym.atan(x_prime2*R*np.sin(theta)/(y_prime2*y4-np.square(x_prime2)*np.cos(theta)))
A=R-z_prime2
B=R-z3
B1=R-z4
def selection(b1,b2):   #select parameters according to different inputs
    if b1!=0:
        global angle,alpha,be,bx
        angle=eval(input("Please input the angle between burgers vector of TED and 11-20 direction"))   #angle between b and 11-20 unit:degree
        alpha=angle*np.pi/180             #angle between b and 11-20 unit:rad
        be,bx=b1*np.cos(alpha),-b1*np.sin(alpha)
        if b2!=0:
            uu=ted1(x,y,z)+tsd(x,y,z)+ted2(x,y,z)
        elif b2==0:
            uu=ted1(x,y,z)+ted2(x,y,z)
    elif b1==0 and b2!=0:
        uu=tsd(x,y,z)
    return uu
def ted1(x,y,z):     #calculation for part of TED component
    global qe,me,lamda,thetaa,k,C
    qe=x_prime2*(1/B1-1/B+(2*z_prime2*np.cos(theta))/np.square(B))
    me=-qe/R-4*(1-v)*x_prime2*np.square(np.cos(theta))/(R*B)
    lamda=(1-2*v)*(sym.log(B1/B))  #lamda in the paper by Shaibani
    thetaa=2*(1-v)*(omega1-omega)
    k=4*(1-v)*(1-2*v)/np.square(np.tan(theta))
    C=8*np.pi*(1-v)/be          #the value of De*2*μ
    uxted=(x_prime2*me+lamda+(2*np.cos(theta)/B)*(z_prime2+2*(1-v)*y3*np.sin(theta))-4*(1-v)*
        np.square(np.sin(theta))+k*(1-np.cos(theta)-np.cos(theta)*sym.log(A)+y_prime2*np.sin(theta)/A+sym.log(B)))/C
    uyted=(y_prime2*me+qe*np.sin(theta)+thetaa*np.cos(theta)+k*(-x_prime2*np.sin(theta)/A+omega*np.cos(theta)))/C
    uzted=(z_prime2*me+qe*np.cos(theta)+thetaa*np.sin(theta)-2*x_prime2*np.cos(theta)*(1/B1+(1-2*v)/B)+k*omega*np.sin(theta))/C
    return np.array([uyted,-uxted,uzted])
def tsd(x,y,z):     #calculation for TSD component
    ms=x_prime2*np.sin(2*theta)/(R*B)
    S=4*np.pi/bs               #S=2*u*s
    uxtsd=(x_prime2*ms+2*y3*np.square(np.cos(theta))/B+2*(1-2*v)*(1/np.tan(theta))*(-1+np.cos(theta)+np.cos(theta)*
          sym.log(A)-y_prime2*np.sin(theta)/A-sym.log(B))-np.sin(2*theta))/S
    uytsd=(y_prime2*ms-2*x_prime2*np.cos(theta)/B-np.sin(theta)*(omega1-omega)+2*(1-2*v)*(1/np.tan(theta))
           *(x_prime2*np.sin(theta)/A-omega*np.cos(theta)))/S
    uztsd=(z_prime2*ms+np.cos(theta)*(omega1-omega)-2*(1-2*v)*omega*np.cos(theta))/S
    return np.array([uytsd,-uxtsd,uztsd])
def ted2(x,y,z):     #calculation for part of TED component
    qx=y4/B1-y3/B-2*z_prime2*y3*np.cos(theta)/np.square(B)
    mx=-qx/R+2*(1-2*v)*y_prime2*np.cos(theta)/(R*B)
    Dx=4*np.pi*(1-v)/bx       #D=D*u
    uxted2=(x_prime2*mx+thetaa+k*(x_prime2*np.tan(theta)/A-omega))/(2*Dx)
    uyted2=(y_prime2*mx+qx*np.sin(theta)-lamda*np.cos(theta)-(2*np.cos(theta)/B)*(z_prime2*np.cos(theta)+(1-2*v)
          *y_prime2*np.sin(theta))+k*(-1+np.cos(theta)-sym.log(A)+y_prime2*np.tan(theta)/A+np.cos(theta)*sym.log(B)))/(2*Dx)
    uzted2=(z_prime2*mx+qx*np.cos(theta)-lamda*np.sin(theta)-2*y4*np.cos(theta)/B1+(4*np.cos(theta)/B)*((1-v)*y_prime2*np.cos(theta)-(1-2*v)*
       z_prime2*np.sin(theta))+k*np.tan(theta)*(np.cos(theta)-sym.log(A)+np.cos(theta)*sym.log(B))+4*(1-v)*np.cos(theta)/np.tan(theta))/(2*Dx)
    return np.array([uyted2,-uxted2,uzted2])
def grad(f,X,Y,Z):    #calculate the diffracted beam
    fx=sym.lambdify([x,y,z],sym.diff(f,x),'numpy')
    fy=sym.lambdify([x,y,z],sym.diff(f,y),'numpy')
    fz=sym.lambdify([x,y,z],sym.diff(f,z),'numpy')
    return np.array([fx(X,Y,Z),fy(X,Y,Z),fz(X,Y,Z)])
nu=np.dot(n,selection(b1,b2)) #dot product of n and u
der=grad(nu,X,Y,Z)
reder=der.reshape(3,-1).transpose().reshape(-1,w2,3)
n3D=np.resize(n,reder.shape)
n1=n3D-reder
n2=np.empty(n1.shape)
for i in range(n1.shape[0]):
    for j in range(n1.shape[1]):
        n2[i,j,:]=n1[i,j,:]/np.linalg.norm(n1[i,j,:],ord=2)
qh0=q0-2*np.dot(q0,n)*n
qh=np.empty(n2.shape)
for i in range(n2.shape[0]):
    for j in range(n2.shape[1]):
        qh[i,j,:]=q0-2*np.dot(q0,n2[i,j,:])*n2[i,j,:]
xc=D/qh0[2]*qh0[0]
yc=D/qh0[2]*qh0[1]
x1=(D/qh[:,:,2]*qh[:,:,0]-xc).reshape(X.shape)+X
y1=(D/qh[:,:,2]*qh[:,:,1]-yc).reshape(Y.shape)+Y
############# ploting
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(10,10),dpi=300)
ax=fig.add_subplot(1,1,1,aspect="equal")
ax.set_xlim(l1,l2)
ax.set_ylim(l3,l4)
point_size = 0.01
pene=np.linspace(0,p,1*w2) [::-1] #penetration depth measured from the crystal surface
for i in range(1*w2):
    color=int(255*(1-np.exp(-u*(1/np.sin(beta-theta0)+1/np.sin(np.pi/2-beta1+theta0))*pene[i]*0.0001)))
    ax.scatter(x1[:,:,i],y1[:,:,i],c=np.full((1,3),color/255),s=point_size,marker='.')
plt.title('Simulated Image of Dislocation')
plt.show()