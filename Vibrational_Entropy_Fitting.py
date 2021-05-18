import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
read output files and rename the columns as "Temperature,Entropy,T*Entropy".
"""
perfect_supercell = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\perfect supercell-entr.txt",
                                sep = "\s+",header = None, names=["Temperature","S","T_S"])
VN = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\VN-entr.txt",
                 sep = "\s+",header = None, names=["Temperature","S","T_S"])
VGa = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\VGa-entr.txt",
                  sep = "\s+",header = None, names=["Temperature","S","T_S"])
Ni = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\Ni-entr.txt",
                 sep = "\s+",header = None, names=["Temperature","S","T_S"])
Gai = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\Gai-entr.txt",
                  sep = "\s+",header = None, names=["Temperature","S","T_S"])
NGa = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\NGa-entr.txt",
                  sep = "\s+",header = None, names=["Temperature","S","T_S"])
GaN = pd.read_csv(r"C:\Users\66407\Desktop\Entropy_data\GaN-entr.txt",
                  sep = "\s+",header = None, names=["Temperature","S","T_S"])

"""
Define a function to return the expression of the fitting curves.
"""
def TdelatS(x,a,b,c,d):
    return a+b*x+c*x**2+d*x**3

"""
read output files and rename the columns as "Temperature,Entropy,T*Entropy"
"""
xdata = perfect_supercell["Temperature"]
VNdata = VN["T_S"]-perfect_supercell["T_S"]
VGadata = VGa["T_S"]-perfect_supercell["T_S"]
Nidata = Ni["T_S"]-perfect_supercell["T_S"]
Gaidata = Gai["T_S"]-perfect_supercell["T_S"]
NGadata = NGa["T_S"]-perfect_supercell["T_S"]
GaNdata = GaN["T_S"]-perfect_supercell["T_S"]

"""
Curve fitting
"""
VNpopt,VNpcov = curve_fit(TdelatS,xdata,VNdata)
VGapopt,VGapcov = curve_fit(TdelatS,xdata,VGadata)
Nipopt,Nipcov = curve_fit(TdelatS,xdata,Nidata)
Gaipopt,Gaipcov = curve_fit(TdelatS,xdata,Gaidata)
NGapopt,NGapcov = curve_fit(TdelatS,xdata,NGadata)
GaNpopt,GaNpcov = curve_fit(TdelatS,xdata,GaNdata)

"""
print the expression of the fitting curves.
"""
print ("N Vacancy: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*VNpopt))
print ("Ga Vacancy: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*VGapopt))
print ("N Interstitial: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*Nipopt))
print ("Ga Interstitial: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*Gaipopt))
print ("N-Ga Substitial: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*NGapopt))
print ("Ga-N Substitial: T*deltaS = {:.3e}+{:.3e}T+{:.3e}T^2+{:.3e}T^3".format(*GaNpopt))

"""
Plot data curve and fitting curve
"""
fig,ax= plt.subplots(figsize = (6,4),dpi = 200)
T = np.linspace(0, 1900, 1900)

plt.plot(xdata,TdelatS(xdata,*VNpopt),ls="--",c="black")
plt.plot(xdata,TdelatS(xdata,*VGapopt),ls="--",c="black")
plt.plot(xdata,TdelatS(xdata,*Nipopt),ls="--",c="black")
plt.plot(xdata,TdelatS(xdata,*Gaipopt),ls="--",c="black")
plt.plot(xdata,TdelatS(xdata,*NGapopt),ls="--",c="black")
plt.plot(xdata,TdelatS(xdata,*GaNpopt),ls="--",c="black")

plt.scatter(T,VN["T_S"]-perfect_supercell["T_S"],label="N Vacancy")
plt.scatter(T,VGa["T_S"]-perfect_supercell["T_S"],label="Ga Vacancy")
plt.scatter(T,Ni["T_S"]-perfect_supercell["T_S"],label="N Interstitial")
plt.scatter(T,Gai["T_S"]-perfect_supercell["T_S"],label="Ga Interstitial")
plt.scatter(T,NGa["T_S"]-perfect_supercell["T_S"],label="N-Ga Substitial")
plt.scatter(T,GaN["T_S"]-perfect_supercell["T_S"],label="Ga-N Substitial")

ax.set_xlim([0,1900])
ax.set_ylim([-3,3])
plt.xlabel("Temperatre(K)",fontsize = 15)
plt.ylabel("$T*deltaS$(eV)",fontsize = 15)
plt.title("Vibrational Entropy")
plt.grid()
plt.legend()
plt.show()

