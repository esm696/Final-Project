# -*- coding: utf-8 -*-
"""
Created on Fri May 14 00:58:39 2021

@author: Dean
"""
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.colors as colors
from scipy.ndimage.filters import gaussian_filter
#---------------------------------------------------------------------------------------------------------------------------
#Input desired file name in text format
File='Vectors.csv'

#Choose desired projection mode. (EQA=Equal Area Projection, SG=Stereographical projection)
Projection_mode='EQA'

#Controls near or far hemisphere projection. (N=Near hemisphere, F=Far hemisphere)
NOF='N'

#Controls the bin size of heatmap, binR/binA controls the radial/angular bin size, respectivly.
binR=300
binA=300

#Controls the projection axis, choose from X, Y and Z
Projection_Axis='Z'

#Optional, plot out the radial scatter plot, slow for large data set (True/False)
scatter=False

#Optional, apply a gaussian filter to the heat map, input the sigma here as well (True/False)
Gaussian_On=True
Sigma=2


#---------------------------------------------------------------------------------------------------------------------------
f1=open(File, "r")
field=[]
c=0
Valid_Plot=True
Valid_Axis=True

if Projection_Axis=='X':
    AXX=1
    AXY=0
    AXZ=0
elif Projection_Axis=='Y':
    AXX=0
    AXY=1
    AXZ=0
elif Projection_Axis=='Z':
    AXX=0
    AXY=0
    AXZ=1  
else:
    Valid_Axis=False
    print('Invalid axis, supported axis are X, Y and Z')



if NOF=='F':
    NF=-1
elif NOF=='N':
    NF=1
else:
    print('Invalid input, running the projection with near hemisphere projection')
    NF=1
        
if Projection_mode=='EQA':
    for i,j in enumerate(f1):
        if i!=0:
            x=float((j.split(',')[4]))
            y=float((j.split(',')[5]))
            z=float((j.split(',')[6]))
            if sqrt(x**2+y**2+z**2)!=0:
                mag=sqrt(x**2+y**2+z**2)
                xu=x/mag
                yu=y/mag
                zu=z/mag                
                cmag=sqrt(xu**2+yu**2+zu**2)
                d=sqrt((AXX*NF-xu)**2+(AXY*NF-yu)**2+(AXZ*NF-zu)**2)
                if Projection_Axis=='X':
                    u1=xu
                    u2=yu
                    u3=zu
                elif Projection_Axis=='Y':
                    u1=yu
                    u2=zu
                    u3=xu 
                elif Projection_Axis=='Z':
                    u1=zu
                    u2=xu
                    u3=yu
                ang=arctan2(u2,u3) 
                if NF==1 and u1<0:
                        continue
                elif NF==-1 and u1>0:
                        continue
                else:
                    field.append([ang,d])
                c=c+1
    print(str(c)+' valid data read')
    
elif Projection_mode=='SG':
    for i,j in enumerate(f1):
        if i!=0:
            x=float((j.split(',')[4]))
            y=float((j.split(',')[5]))
            z=float((j.split(',')[6]))
            if sqrt(x**2+y**2+z**2)!=0:
                mag=sqrt(x**2+y**2+z**2)
                xu=x/mag
                yu=y/mag
                zu=z/mag
                cmag=sqrt(xu**2+yu**2+zu**2)
                if Projection_Axis=='X':
                    u1=xu
                    u2=yu
                    u3=zu
                elif Projection_Axis=='Y':
                    u1=yu
                    u2=zu
                    u3=xu 
                elif Projection_Axis=='Z':
                    u1=zu
                    u2=xu
                    u3=yu
                AN=sqrt((xu-(-AXX*NF))**2+(yu-(-AXY*NF))**2+(zu-(-AXZ*NF))**2)
                if Projection_Axis=='X':
                    AZ=abs(xu-(-AXX*NF))   
                if Projection_Axis=='Y':
                    AZ=abs(yu-(-AXY*NF))  
                if Projection_Axis=='Z':
                    AZ=abs(zu-(-AXZ*NF))         

                if AN!=0:
                    if NF==1 and u1<0:
                            continue
                    elif NF==-1 and u1>0:
                            continue
                    else:
                        cosA=AZ/AN
                        pang=arccos(cosA)
                        azu=arctan2(u2,u3)  
                        ANN=2/cos(pang)
                        R=sqrt(ANN**2-4)        
                        field.append([azu,R])
                        c=c+1
    print(str(c)+' valid data read')
else:
    Valid_Plot=False
    print('Not supported prjection method')



if scatter==True and Valid_Plot==True:
    print('\n'+'Preparing polar plot...')
    fig0 = plt.figure()    
    for i,j in enumerate(field):
        if (i+1)%1000==0:
            if i!=0:
                print(str(i+1)+' data processed')
        elif i+1==(c):
            print('All('+str(i+1)+ ') data processed')
            print('Processing polar plot')      
        plt.polar((j[0]+pi),j[1], 'g.')

            
    ax0 = fig0.add_subplot(111, projection='polar')
    if Projection_mode=='EQA':
        ax0.set_ylim(0,sqrt(2))
    else:
        ax0.set_ylim(0,2)
    ax0.set_yticklabels([])
    plt.show()
    
if Valid_Plot==True:
    fig = plt.figure()
    rad=[]
    deg=[]
    for i,j in enumerate(field) :
        rad.append(j[1])
        deg.append(j[0])

    rad2=array(rad)    
    deg2=array(deg)

    rbins = linspace(0,rad2.max(), binR)
    abins = linspace(0,2*pi, binA)
    hist, _, _ = histogram2d(deg2, rad2, bins=(binA,binR),density=True)
    if Gaussian_On==True:
        hist=gaussian_filter(hist, sigma=Sigma)
    A, R = meshgrid(abins, rbins)
    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])   
    pc = ax.pcolormesh(A, R, hist.T,cmap="jet",norm=colors.Normalize(vmin=0, vmax=0.5, clip=True))  