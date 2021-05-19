import collections
from itertools import islice
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib import rc
"""
This file is to draw a bar plot for the Voronoi statistics 
for 4 glass particles, i.e 4 bar IGC particle, 6 bar IGC particle,
10 bar IGC particle and melting-quenching bulk metallic glass(reference)
"""

rc('font',**{'family':'serif','serif':['Arial']}) #set the font


def take(n,iterable):

    return list(islice(iterable,n)) # take elements from a list("iterable") from 0 to n
#read files

with open("voronoi_MQ","r") as file1:
    voronoi1 = file1.readlines();
file1.close()

with open("voronoi_IGC4","r") as file2:
    voronoi2 = file2.readlines();
file2.close()

with open("voronoi_IGC6","r") as file3:
    voronoi3 = file3.readlines();
file3.close()

with open("voronoi_IGC10","r") as file4:
    voronoi4 = file4.readlines();
file4.close()

#count numbers and sort
vorstatics1 = dict() #using a dictinary to store unique Voronoi index and the frquency
vorstatics2 = dict()
vorstatics3 = dict()
vorstatics4 = dict()

for value in voronoi1:
    if not vorstatics1.keys().__contains__(value):
        vorstatics1.update({value:voronoi1.count(value)})
#vorstatics1 = dict()
for value in voronoi2:
    if not vorstatics2.keys().__contains__(value):
        vorstatics2.update({value:voronoi2.count(value)})
for value in voronoi3:
    if not vorstatics3.keys().__contains__(value):
        vorstatics3.update({value:voronoi3.count(value)})
for value in voronoi4:
    if not vorstatics4.keys().__contains__(value):
        vorstatics4.update({value:voronoi4.count(value)})

#sort from most prevelent to the least
sort_vori1 = sorted(vorstatics1.items(),key = lambda x:x[1],reverse=True)
sort_vori2 = sorted(vorstatics2.items(),key = lambda x:x[1],reverse=True)
sort_vori3 = sorted(vorstatics3.items(),key = lambda x:x[1],reverse=True)
sort_vori4 = sorted(vorstatics4.items(),key = lambda x:x[1],reverse=True)
#Take the first 10
vor_top101 = take(10,sort_vori1)
vor_top102 = take(10,sort_vori2)
vor_top103 = take(10,sort_vori3)
vor_top104 = take(10,sort_vori4)

#vor_top10_paper = [vor_top10[0],vor_top10[1],vor_top10[6],vor_top10[2],vor_top10[3],\
#                   vor_top10[8],vor_top10[5],vor_top10[4],vor_top10[9],vor_top10[7]]

#using sorted vor_top10
#countnum = np.array([value[1] for value in vor_top10])
#fraction = (countnum/float(len(data)))*100
#keys = [value[0] for value in vor_top10]

#using the order in paper
countnum1 = np.array([value[1] for value in vor_top101])
fraction1 = np.round((countnum1/float(len(voronoi1)))*100,5)#calculate the concentration of Voronoi index
keys1 = [value[0] for value in vor_top101]
print(keys1)
#sort
#countnum2 = np.array([value[1] for value in vor_top102])
#fraction2 = (countnum2/float(len(voronoi2)))*100
#keys2 = [value[0] for value in vor_top102]
a = 0
keys2 = []
fraction2 = []
print("2")
for i in range(0,len(sort_vori2)):
    for j in range(0,len(vor_top101)):
        if sort_vori2[i][0]==vor_top101[j][0]:
            fraction2.append((round((float(sort_vori2[i][1]) / float(len(voronoi2))*100), 5)))  #check: let round to 5 digits
            keys2.append(sort_vori2[i][0])
            a = a + 1

    if a==len(vor_top101):
        break;

for i in range(0,len(keys2)):
    if keys2[i] != keys1[i]:
        t = keys2[i]
        t_value = fraction2[i]
        num = keys2.index(keys1[i])
        keys2[i]=keys2[num]
        fraction2[i]=fraction2[num]
        keys2[num]=t
        fraction2[num]=t_value

a = 0
keys3 = []
fraction3 = []
print(3)
for i in range(0,len(sort_vori3)):
    for j in range(0,len(vor_top101)):
        if sort_vori3[i][0]==vor_top101[j][0]:
            fraction3.append((round((float(sort_vori3[i][1]) / float(len(voronoi3))*100), 5)))   #check: let round to 5 digits
            keys3.append(sort_vori3[i][0])
            a = a + 1

    if a==len(vor_top101):
        break;

for i in range(0,len(keys3)):
    if keys3[i] != keys1[i]:
        t = keys3[i]
        t_value = fraction3[i]
        num = keys3.index(keys1[i])
        keys3[i] = keys3[num]
        fraction3[i] = fraction3[num]
        keys3[num] = t
        fraction3[num] = t_value

print(4)
a = 0
keys4 = []
fraction4 = []
for i in range(0,len(sort_vori4)):
    for j in range(0,len(vor_top101)):
        if sort_vori4[i][0]==vor_top101[j][0]:
            fraction4 .append(round((float(sort_vori4[i][1]) / float(len(voronoi4)))*100,5))   #check: let round to 5 digits
            keys4.append(sort_vori4[i][0])
            a = a + 1
    if a==len(vor_top101):
        break;

for i in range(0,len(keys4)):
    if keys4[i] != keys1[i]:
        t = keys4[i]
        t_value = fraction4[i]
        num = keys4.index(keys1[i])
        keys4[i] = keys4[num]
        fraction4[i] = fraction4[num]
        keys4[num] = t
        fraction4[num] = t_value

print(keys1)
print(keys2)
print(keys3)
print(keys4)


#plot the bar chart
width = 0.2 #set width
x1=np.arange(len(keys1))-width*1.5 #set the center of bar
x2=np.arange(len(keys2))-width/2
x3=np.arange(len(keys3))+width/2
x4=np.arange(len(keys4))+width*1.5

#x2=np.arange(len(keys2))-width
#x3=np.arange(len(keys3))
#x4=np.arange(len(keys4))+width


fig,ax = plt.subplots()
ax.tick_params(which='both',right=True, top=True)

ax.bar(x1,fraction1,width,align='center', facecolor='black',alpha=0.8,zorder=3,label='melting-quenching particle')
ax.bar(x2,fraction2,width,align='center', facecolor='green',alpha=0.8,zorder=3,label='IGC particle at 4 bars')
ax.bar(x3,fraction3,width,align='center', facecolor='blue',alpha=0.8,zorder=3,label='IGC particle at 6 bars')
ax.bar(x4,fraction4,width,align='center', facecolor='orange',alpha=0.8,zorder=3,label='IGC particle at 10 bars')

plt.xticks(x2+width*0.5,keys1,rotation=45)
plt.ylim([0,20])
plt.ylabel('Fraction of atoms(%)',fontsize=20)
plt.xlabel('Polyhedron Index',fontsize=20)
l1 = ax.legend(loc='upper right', frameon=False, fontsize=16)
ax.yaxis.set_major_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.tick_params(which='major',direction='in',labelsize=16,length=4,width=2,color=[0,0,0],zorder=100)
ax.tick_params(which='minor',direction='in',length=3,width=2,color=[0,0,0],zorder=100)
for side in ax.spines.keys():
    ax.spines[side].set_linewidth(2)

plt.tight_layout()
plt.show()


