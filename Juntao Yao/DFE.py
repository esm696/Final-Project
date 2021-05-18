import numpy as np
import matplotlib.pyplot as plt

"""
Create a data_list describing six different defects' data: "label","transition_level","charge","DFE".
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
Take the N vacancy as an example.
"""
def N_rich(axes):
    for element in data_list_N:
        def f(x):
            len_ = len(element["transition_level"])
            for i in range(len_):
                if element["label"] == "VN":
                    if x <= element["transition_level"][i]:
                        """
                        Assuming that x=0.5, then 0.5<=0.68 (The first object in charge array,i=0).
                        So, DFE_plot will return ["charge"][3-1-0] * x + ["DFE"][3-1-0] - 8.29
                                                 ["charge"][2]=3
                                                 
                        Assuming that x=0.7, then 0.7>0.68 and 0.7<=3.17 (The second object in charge array,i=1).
                        So, DFE_plot will return ["charge"][3-1-1] * x + ["DFE"][3-1-1] - 8.29
                                                 ["charge"][1]=1
                                                 
                        Assuming that x=3.20, then 3.20>3.17 and 3.20<=3.50 (The third object in charge array,i=2).
                        So, DFE_plot will return ["charge"][3-1-2] * x + ["DFE"][3-1-2] - 8.29
                                                 ["charge"][0]=0
                                                 
                        8.29 is the chemical potential for one N atom in N-rich condition.
                        """
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 8.29
                        return DFE_plot
                elif element["label"] == "Ni":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 8.29
                        return DFE_plot
                elif element["label"] == "VGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 3.84
                        """
                        3.84 is the chemical potential for one Ga atom in N-rich condition.
                        """
                        return DFE_plot
                elif element["label"] == "Gai":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 3.84
                        return DFE_plot
                elif element["label"] == "NGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 3.84 + 8.29
                        return DFE_plot
                elif element["label"] == "GaN":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 3.84 - 8.29
                        return DFE_plot

        axes[0].plot(bandgap, np.vectorize(f)(bandgap), label=element["label"])

def N_poor(axes):
    for element in data_list_N:
        def f(x):
            len_ = len(element["transition_level"])
            for i in range(len_):
                if element["label"] == "VN":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 9.21
                        """
                        9.21 is the chemical potential for one N atom in N-poor condition.
                        """
                        return DFE_plot
                elif element["label"] == "Ni":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 9.21
                        return DFE_plot
                elif element["label"] == "VGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 2.92
                        """
                        2.92 is the chemical potential for one Ga atom in N-poor condition.
                        """
                        return DFE_plot
                elif element["label"] == "Gai":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 2.92
                        return DFE_plot
                elif element["label"] == "NGa":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] - 2.92 + 9.21
                        return DFE_plot
                elif element["label"] == "GaN":
                    if x <= element["transition_level"][i]:
                        DFE_plot = element["charge"][len_ - 1 - i] * x + element["DFE"][len_ - 1 - i] + 2.92 - 9.21
                        return DFE_plot

        axes[1].plot(bandgap, np.vectorize(f)(bandgap), label=element["label"])


bandgap = np.linspace(0, 3.5, 100)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 5), dpi=200)

N_rich(axes)
N_poor(axes)

axes[0].set_yticks(np.arange(-2, 11))
axes[1].set_yticks(np.arange(-2, 11))
axes[0].set_xlim([0, 3.5])
axes[1].set_xlim([0, 3.5])
axes[0].set_xlabel("Fermi level(eV)")
axes[1].set_xlabel("Fermi level(eV)")
axes[0].set_ylabel("Defects Formation Energy(eV)")
axes[0].set_title("N-rich")
axes[1].set_title("N-poor")
axes[0].legend()
axes[1].legend()
axes[0].grid()
axes[1].grid()
plt.show()
