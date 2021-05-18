import numpy as np
import matplotlib.pyplot as plt
from Vibrational_Entropy_Fitting import VNpopt,VGapopt,Nipopt,Gaipopt,NGapopt,GaNpopt,TdelatS
"""
Combine DFT and thermodynamic results, adding the contribution of the vibrational entropy.
"""
data_list_N = [
    {
        "label": "VN",
        "transition_level": np.array([0.68, 3.17, 3.5]),
        "charge": np.array([0, 1, 3]),
        "DFE": np.array([11.32, 8.15, 6.78]),
    },
    {
        "label": "Ni",
        "transition_level": np.array([0.77, 1.80, 2.20, 3.28, 3.50]),
        "charge": np.array([-1, 0, 1, 2, 3]),
        "DFE": np.array([0.74, -2.54, -4.73, -6.53, -7.30]),
    },
    {
        "label": "VGa",
        "transition_level": np.array([1.88, 2.10, 3.13, 3.50]),
        "charge": np.array([-3, -2, -1, 0]),
        "DFE": np.array([17.54, 14.4, 12.3, 10.42, 9.79]),
    },
    {
        "label": "Gai",
        "transition_level": np.array([2.43, 2.83, 3.50]),
        "charge": np.array([1, 2, 3]),
        "DFE": np.array([2.44, -0.39, -2.82]),
    },
    {
        "label": "NGa",
        "transition_level": np.array([1.70, 2.67, 3.50]),
        "charge": np.array([0, 1, 2]),
        "DFE": np.array([2.59, -0.09, -1.79]),
    },
    {
        "label": "GaN",
        "transition_level": np.array([1.56, 1.60, 2.13, 2.50, 3.23, 3.50]),
        "charge": np.array([-1, 0, 1, 2, 3, 4]),
        "DFE": np.array([15.53, 12.29, 9.79, 7.66, 6.11, 4.5]),
    }
]
"""
Set a temperature
"""
T = 1800
def N_rich(axes):
    for element in data_list_N:
        def f(x):
            len_ = len(element["transition_level"])
            for i in range(len_):
                if element["label"] == "VN":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 8.29
                        return DFE_plot + TdelatS(T,*VNpopt)
                    """
                    Call TdeltaST to adding T*deltaS
                    """
                elif element["label"] == "Ni":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 8.29
                        return DFE_plot + TdelatS(T,*Nipopt)
                elif element["label"] == "VGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 3.84
                        return DFE_plot + TdelatS(T,*VGapopt)
                elif element["label"] == "Gai":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 3.84
                        return DFE_plot + TdelatS(T,*Gaipopt)
                elif element["label"] == "NGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 3.84 + 8.29
                        return DFE_plot + TdelatS(T,*NGapopt)
                elif element["label"] == "GaN":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 3.84 - 8.29
                        return DFE_plot + TdelatS(T,*GaNpopt)

        axes.plot(bandgap, np.vectorize(f)(bandgap), label=element["label"])

bandgap = np.linspace(0, 3.5, 100)
fig, axes = plt.subplots(figsize=(2.5,5), dpi=200)
N_rich(axes)
axes.set_yticks(np.arange(-2, 11))
axes.set_xlim([0, 3.5])
axes.set_xlabel("Fermi level(eV)")
axes.set_ylabel("Defects Formation Energy(eV)")
axes.legend()
axes.grid()
plt.show()