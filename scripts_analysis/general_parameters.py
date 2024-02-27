import sys,os
import pickle as pkl
#from scripts.parse_tools import main_parser,collect_results
#from scripts.plot_tools import *
from general_tools import lin_fun,quad_fun
from ase.io import read, write
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
import matplotlib
from matplotlib import gridspec

matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['figure.figsize'] = (8,6)
markersize=10
from matplotlib.ticker import FormatStrFormatter


home=os.getcwd()
SHE_potential=4.4
Aperatm=6.928203230275509
refs={'H': -8.009969/2.,'Na':-21.349228,'K': -40.708912,'Cs':-183.496667,'Li':-5.595442,'Mg':-47.498754/2.}
# Redox potentials come from Vanýsek, Petr (2011). "Electrochemical Series".
# In Haynes, William M. (ed.). CRC Handbook of Chemistry and Physics (92nd ed.). CRC Press. pp. 5–80–9. ISBN 978-1-4398-5512-6.
redoxpots={'Na':SHE_potential-2.71,'K':SHE_potential-2.931,'H':SHE_potential,'Cs':SHE_potential-3.026, 'Li': SHE_potential-3.04,'Mg':SHE_potential-2.362}
elneg={'H':2.2,'Li':0.98,'Na':0.93,'K':0.82,'Cs':0.79,'Mg':1.31}
sizes=['2x2','sqrt7xsqrt7','3x3','sqrt13xsqrt13','4x4']


catcolors=['k','r','b','g','y','brown','orange','r','r']
colors=[cm.nipy_spectral(i/float(len(sizes)+1)) for i in range(len(sizes)+1)]
del colors[0]
markers=['o','X','d','^','>','<']


coverages_for_legends={'sqrt3xsqrt3':'1/3',
            '2x2': r'$\frac{1}{4}$',
            'sqrt7xsqrt7': r'$\frac{1}{7}$',
            '3x3': r'$\frac{1}{9}$',
            'sqrt13xsqrt13': r'$\frac{1}{13}$',
            '4x4': r'$\frac{1}{16}$'}
coverages={'sqrt3xsqrt3':0.33,
            '2x2': 0.25,
            'sqrt7xsqrt7': 0.14,
            '3x3': 0.11,
            'sqrt13xsqrt13': 0.077,
            '4x4': 0.0625}
