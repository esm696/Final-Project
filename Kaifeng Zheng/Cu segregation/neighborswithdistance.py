import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline

#read files
with open('10barNND','r') as file1:
    data1 = file1.readlines()

with open('6barNND','r') as file2:
    data2 = file2.readlines()

with open('4barNND','r') as file3:
    data3 = file3.readlines()


id10 = []    # Particle identifier
coord10 = [] # number of neighbors with in a spherical region with 10 angstrom radius 
type10 = []  # 3: Zr 2: Cu


id6 = []# Particle identifier
coord6 = []# number of neighbors with in a spherical region with 10 angstrom radius 
distance6 = []# center to atom distance
type6 = [] # 3: Zr 2: Cu

id4 = []# Particle identifier
coord4 = []# number of neighbors with in a spherical region with 10 angstrom radius 
distance4 = []# center to atom distance
type4 = [] # 3: Zr 2: Cu

for i in range(1,len(data1)):
    temp = data1[i].split()
    density10 = (float(temp[1])/((4.0/3.0)*np.pi*10**3))*100 # calculate number density in the region with 10 angstrom radius
    id10.append(int(temp[0]))# get particle identifier
    coord10.append(density10)# get nnumber of neighbors
    type10.append(int(temp[2]))# get type of atom

for i in range(1,len(data2)):
    temp = data2[i].split()
    density6 = (float(temp[1]) / ((4.0 / 3.0) * np.pi*10**3))*100 # calculate number density in the region with 10 angstrom radius
    id6.append(int(temp[0]))# get particle identifier
    coord6.append(density6)# get nnumber of neighbors
    distance6.append(float(temp[2]))# get center to atom distance
    type6.append(int(temp[3]))# get type 

for i in range(1,len(data3)):
    temp = data3[i].split()
    density4 = (float(temp[1]) / ((4.0 / 3.0) * np.pi*10**3))*100 # calculate number density in the region with 10 angstrom radius
    id4.append(int(temp[0]))# get particle identifier
    coord4.append(density4)# get nnumber of neighbors
    distance4.append(float(temp[2]))# get center to atom distance
    type4.append(int(temp[3]))# get type 


mincoo10 = min(coord10)#find min number of neighbors
maxcoo10 = max(coord10)#find max number of neighbors


co2dis6 = dict() #{neighbors:distance}
mincoo6 = min(coord6)
maxcoo6 = max(coord6)

co2dis4 = dict()
mincoo4 = min(coord4)
maxcoo4 = max(coord4)

# number of neighbors vs atom-center distance
for j in range(0,len(id6)):
    if co2dis6.__contains__(coord6[j]):
        co2dis6[coord6[j]].append(distance6[j])
    else:
        co2dis6[coord6[j]]=list()
        co2dis6[coord6[j]].append(distance6[j])
#sort the dict 
keys6 = sorted(co2dis6)

co2dis_sort6 = dict()
for i in keys6[5:149]:
    co2dis_sort6[i]=round(sum(co2dis6[i])/len(co2dis6[i]),2) # get the everage distance


for j in range(0,len(id4)):
    if co2dis4.__contains__(coord4[j]):
        co2dis4[coord4[j]].append(distance4[j])
    else:
        co2dis4[coord4[j]]=list()
        co2dis4[coord4[j]].append(distance4[j])
keys4 = sorted(co2dis4)

co2dis_sort4 = dict()
for i in keys4[5:149]:
    co2dis_sort4[i]=round(sum(co2dis4[i])/len(co2dis4[i]),2)

# fitting the data using linear fitting
fitting6 = np.polyfit(list(co2dis_sort6.keys()),list(co2dis_sort6.values()),1)
p6 = np.poly1d(fitting6)
fitting4 = np.polyfit(list(co2dis_sort4.keys()),list(co2dis_sort4.values()),1)
p4 = np.poly1d(fitting4)

# plot
fig, ax1 = plt.subplots()

ax1.plot(distance6,coord6,'.',color=[0,1,0],alpha=0.2,label='6 bars')
ax1.plot(distance4,coord4,'b.',alpha=0.1,label='4 bars')
ax1.plot([0,14],[6.37207,6.37207],':',color=[1,0,0.498],linewidth=2)
ax1.text(0.1,6.43,r'$\rho_{bulk}=6.37208\times10^{-2}\AA^{-3}$',color='r',fontsize=18)
ax1.plot(p6(list(co2dis_sort6.keys())),list(co2dis_sort6.keys()),'-',color='red',linewidth=2,label='trend for 6 bars')
ax1.plot(p4(list(co2dis_sort4.keys())),list(co2dis_sort4.keys()),'-',color='orange',linewidth=2,label='trend for 4 bars')
#ax1.plot(distance4,p4(distance4),'--',color='black')

l1 = ax1.legend(loc='lower right', frameon=False, fontsize=18)
for lh in l1.legendHandles:
    lh._legmarker.set_alpha(1)
ax1.tick_params(which='both',right=True, top=True)
ax1.xaxis.set_major_locator(plt.MultipleLocator(2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax1.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
#bp2 = ax3.boxplot(coo4,sym = 'k+', positions=pos,widths=0.25)
ax1.tick_params(which='major',direction='in',labelsize=16,length=4,width=2,color=[0,0,0])
ax1.tick_params(which='minor',direction='in',labelsize=16,length=3,width=2,color=[0,0,0])
plt.xlabel('Distance from particle surface ($\AA$)',fontsize=20)
plt.ylabel(r'Local density ($10^{-2}\AA^{-3}$)',fontsize=20)

for side in ax1.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax1.spines[side].set_linewidth(2)
plt.xlim(0,14)



## Cu segregation test

co2element10 = dict() #{number of neighbors: type of atoms}
co2element6 = dict()
co2element4 = dict()
for j in range(0,len(id10)):
    if co2element10.__contains__(coord10[j]):
        co2element10[coord10[j]].append(type10[j])
    else:
        co2element10[coord10[j]]=list()
        co2element10[coord10[j]].append(type10[j])
keyselement10 = sorted(co2element10) #sort using number of neighbors


for j in range(0,len(id6)):
    if co2element6.__contains__(coord6[j]):
        co2element6[coord6[j]].append(type6[j])
    else:
        co2element6[coord6[j]]=list()
        co2element6[coord6[j]].append(type6[j])
keyselement6 = sorted(co2element6)

for j in range(0,len(id4)):
    if co2element4.__contains__(coord4[j]):
        co2element4[coord4[j]].append(type4[j])
    else:
        co2element4[coord4[j]]=list()
        co2element4[coord4[j]].append(type4[j])
keyselement4 = sorted(co2element4)

################################################TRANFORMATION########################################################################
element10 = dict()
element6 = dict()
element4 = dict()
T_k10 = []
T_k6 = []
T_k4 = []
y10 = []
y6 = []
y4 = []
for i in keyselement10:
    T_k10.append(p6(i)) #add the linear fitting parameters
for i in keyselement6:
    T_k6.append(p6(i))
for i in keyselement4:
    T_k4.append(p6(i))

##############################################END TRANSFORMATION######################################################################
for i in keyselement10:
    num = 0
    for j in range(0,len(co2element10[i])):
        if co2element10[i][j]==2:
            num = num + 1
    element10[i]=round((num/len(co2element10[i])),2)

for i in keyselement6:
    num = 0
    for j in range(0,len(co2element6[i])):
        if co2element6[i][j]==2:
            num = num + 1
    element6[i]=round((num/len(co2element6[i])),2)

for i in keyselement4:
    num = 0
    for j in range(0,len(co2element4[i])):
        if co2element4[i][j]==2:
            num = num + 1
    element4[i]=round((num/len(co2element4[i])),2)
"""
This part is to claculate average position
"""

avex10 = []
avey10 = []
avex6 = []
avey6 = []
avex4 = []
avey4 = []
for i in range(0,len(T_k10),3):
    sumx = 0
    sumy = 0
    if i+3>len(T_k10):
        break
    for j in range(0,3):
        sumx = sumx+T_k10[i+j]
        sumy = sumy+list(element10.values())[i+j]
    avex10.append(round(sumx/3,10))
    avey10.append(round(sumy/3,10))

for i in range(0,len(T_k6),3):
    sumx = 0
    sumy = 0
    if i+3>len(T_k6):
        break
    for j in range(0,3):
        sumx = sumx+T_k6[i+j]
        sumy = sumy+list(element6.values())[i+j]
    avex6.append(round(sumx/3,10))
    avey6.append(round(sumy/3,10))

for i in range(0,len(T_k4),3):
    sumx = 0
    sumy = 0
    if i+3>len(T_k4):
        break
    for j in range(0,3):
        sumx = sumx+T_k4[i+j]
        sumy = sumy+list(element4.values())[i+j]
    avex4.append(round(sumx/3,10))
    avey4.append(round(sumy/3,10))

sp1 = UnivariateSpline(avex10[7:],avey10[7:],k=5,s=900)
sp2 = UnivariateSpline(avex6[7:],avey6[7:],k=5,s=700)
sp3 = UnivariateSpline(avex4[7:],avey4[7:],k=5,s=900)
fig, ax2 = plt.subplots()


ax2.tick_params(which='both',right=True, top=True)
ax2.plot(T_k10,list(element10.values()),'.',color='cyan',label = '10 bars')
ax2.plot(T_k6,list(element6.values()),'.',color = [0,1,0],label = '6 bars')
ax2.plot(T_k4,list(element4.values()),'.',color = 'blue', label = '4 bars')
ax2.plot(avex10,sp1(avex10),'-',color =[1,0,1],linewidth=2,label = 'trend for 10 bars')
ax2.plot(avex6,sp2(avex6),'-',color = 'red',linewidth=2,label = 'trend for 6 bars')
ax2.plot(avex4,sp3(avex4),'-',color ='orange',linewidth=2,label = 'trend for 4 bars')

ax2.plot([0,10],[0.64,0.64],':',color=[1,0,0.5],linewidth=2)
plt.xlim(0,8)
plt.ylim(0,1)

l1 = ax2.legend(loc='lower left', frameon=False, fontsize=18)
for lh in l1.legendHandles:
    lh._legmarker.set_alpha(1)

ax2.xaxis.set_major_locator(plt.MultipleLocator(1))
ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax2.yaxis.set_major_locator(plt.MultipleLocator(0.16))
ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.08))

ax2.tick_params(which='major',direction='in',labelsize=16,length=4,width=2,color=[0,0,0])
ax2.tick_params(which='minor',direction='in',labelsize=16,length=3,width=2,color=[0,0,0])
plt.xlabel('Distance from particle surface ($\AA$)',fontsize=20)
plt.ylabel('Cu atomic fraction',fontsize=20)

for side in ax2.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax2.spines[side].set_linewidth(2)



fig,ax3 = plt.subplots()

avex10 = []
avey10 = []
avex6 = []
avey6 = []
avex4 = []
avey4 = []

for i in range(0,len(list(element10.keys())),3):
    sumx = 0
    sumy = 0
    if i+3>len(list(element10.keys())):
        break
    for j in range(0,3):
        sumx = sumx+list(element10.keys())[i+j]
        sumy = sumy+list(element10.values())[i+j]
    avex10.append(round(sumx/3,10))
    avey10.append(round(sumy/3,10))

for i in range(0,len(list(element6.keys())),3):
    sumx = 0
    sumy = 0
    if i+3>len(list(element6.keys())):
        break
    for j in range(0,3):
        sumx = sumx+list(element6.keys())[i+j]
        sumy = sumy+list(element6.values())[i+j]
    avex6.append(round(sumx/3,10))
    avey6.append(round(sumy/3,10))

for i in range(0,len(list(element4.keys())),3):
    sumx = 0
    sumy = 0
    if i+3>len(list(element4.keys())):
        break
    for j in range(0,3):
        sumx = sumx+list(element4.keys())[i+j]
        sumy = sumy+list(element4.values())[i+j]
    avex4.append(round(sumx/3,10))
    avey4.append(round(sumy/3,10))



sp1 = UnivariateSpline(avex10[7:],avey10[7:],k=5,s=800)
sp2 = UnivariateSpline(avex6[7:],avey6[7:],k=5,s=600)
sp3 = UnivariateSpline(avex4[7:],avey4[7:],k=5,s=800)

ax3.tick_params(which='both',right=True, top=True)
ax3.plot(list(element10.keys()),list(element10.values()),'.',color = 'cyan',label = '10 bars')
ax3.plot(list(element6.keys()),list(element6.values()),'.',color = [0,1,0],label = '6 bars')
ax3.plot(list(element4.keys()),list(element4.values()),'.',color = 'blue', label = '4 bars')

ax3.plot(avex10,sp1(avex10),'-',color =[1,0,1],linewidth=2,label = 'trend for 10 bars')
ax3.plot(avex6,sp2(avex6),'-',color = 'red',linewidth=2,label = 'trend for 6 bars')
ax3.plot(avex4,sp3(avex4),'-',color ='orange',linewidth=2,label = 'trend for 4 bars')


ax3.plot([0,10],[0.64,0.64],':',color=[1,0,0.5],linewidth=2)
plt.ylim(0,1)
plt.xlim(2.1,6)
l1 = ax3.legend(loc='lower left', frameon=False, fontsize=18)
for lh in l1.legendHandles:
    lh._legmarker.set_alpha(1)

ax3.xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax3.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax3.yaxis.set_major_locator(plt.MultipleLocator(0.16))
ax3.yaxis.set_minor_locator(plt.MultipleLocator(0.08))

ax3.tick_params(which='major',direction='in',labelsize=16,length=4,width=2,color=[0,0,0])
ax3.tick_params(which='minor',direction='in',labelsize=16,length=3,width=2,color=[0,0,0])

plt.xlabel(r'Local density ($10^{-2}\AA^{-3}$)',fontsize=20)
plt.ylabel('Cu atomic fraction',fontsize=20)

for side in ax3.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax3.spines[side].set_linewidth(2)
