from ovito.data import *
import numpy as np
from itertools import groupby
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from os import path
import os
"""
This is writing in OVITO, please put the script in "Python script" modifier
"""
def modify(frame, data):
    #plt.clf()
    print("The input contains %i particles." % data.particles.count)
    atomN = data.particles.count
    atomLayer = []
    voronoi=data.particles['Voronoi Index'].array
    print(len(voronoi))
    ParticleTy = data.particles['Particle Type'].array
    with open("D:\\Voronoi_text",'w') as file1: #change the path
        for i in range(0,len(voronoi)):
            if ParticleTy[i] == 2:
                file1.write('Cu<'+str(voronoi[i][2])+' '+str(voronoi[i][3])+' '+str(voronoi[i][4])+' '+str(voronoi[i][5])+'>\n')
            if ParticleTy[i] == 3:
                file1.write('Zr<'+str(voronoi[i][2])+' '+str(voronoi[i][3])+' '+str(voronoi[i][4])+' '+str(voronoi[i][5])+'>\n')
    file1.close()
    print(ParticleTy)