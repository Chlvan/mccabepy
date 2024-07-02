#%%
import matplotlib.pyplot as plt
import numpy as np
import scipy
import cython
import pandas as pd
from phasepy import mixture, component, preos
from phasepy.equilibrium import bubblePy, bubbleTy, dewPx, dewTx

# Reading the database
df = pd.read_excel('All Critical Properties and Acentric Factors.xlsx')

# Extract Tc, Pc, and omega from the database
Tcmat = df['Critical Temperature, Tc (K)']
Pcmat = df['Critical Pressure, Pc (Pa)']
omegamat = df['Acentric Factor, ω (unitless)']
Tcindex = pd.Index(Tcmat)
Pcindex = pd.Index(Pcmat)
omegaindex = pd.Index(omegamat)

#extract the materials or substances, CAS name, iupac name, and SMILES, and CAS registry number columns
materialmat = df['material or substance name']
CASmat = df['CAS name']
iupacmat = df['IUPAC name']
SMILESmat = df['SMILES']
CASnummat = df['CAS Registry No.']

# Lowercase all of the strings in the columns to make inputs case insensitive
materialmat = materialmat.str.lower()
CASmat = CASmat.str.lower()
iupacmat = iupacmat.str.lower()
SMILESmat = SMILESmat.str.lower()
CASnummat = CASnummat.str.lower()

# Index all of the columns
materialindex = pd.Index(materialmat)
CASindex = pd.Index(CASmat)
iupacindex = pd.Index(iupacmat)
SMILESindex = pd.Index(SMILESmat)
CASnumindex = pd.Index(CASnummat)

# Converting Pa to bar for the Pci and Pcj values
Pcmat = Pcmat/10**5 # bar

def VLEPR(compi, compj, T = 298, P = 1): # T in K and P in bar
    # Inputs, later turned into parameters in a function
    # P = 101325 # Pa, guess for Pxy and P for Txy and xy
    # T = 315 # K, guess for Txy and T for Pxy and xy
    # compi = 'methanol' # doesn't have to be iupac
    # compj = 'water' # doesn't have to be iupac

    # Make all components lowercase
    compi = compi.lower()
    compj = compj.lower()

    # Find the index of the components
    if compi in iupacindex:
        i = iupacindex.get_loc(compi)
    elif compi in materialindex:
        i = materialindex.get_loc(compi)
    elif compi in CASindex:
        i = CASindex.get_loc(compi)
    elif compi in SMILESindex:
        i = SMILESindex.get_loc(compi)
    elif compi in CASnumindex:
        i = CASnumindex.get_loc(compi)
    else:
        print('First component not found in database')
    if compj in iupacindex:
        j = iupacindex.get_loc(compj)
    elif compj in materialindex:
        j = materialindex.get_loc(compj)
    elif compj in CASindex:
        j = CASindex.get_loc(compj)
    elif compj in SMILESindex:
        j = SMILESindex.get_loc(compj)
    elif compj in CASnumindex:
        j = CASnumindex.get_loc(compj)
    else:
        print('Second component not found in database')

    # Extract the Tc, Pc, and omega values for the components
    Pci = Pcmat[i] # Pa
    Tci = Tcmat[i] # K
    omegai = omegamat[i] # unitless
    Pcj = Pcmat[j] # Pa
    Tcj = Tcmat[j] # K
    omegaj = omegamat[j] # unitless

    # Create the components
    compi = component(name=iupacmat[i], Tc=Tci, Pc=Pci, w=omegai)
    compj = component(name=iupacmat[j], Tc=Tcj, Pc=Pcj, w=omegaj)
    mix = mixture(compi, compj) # mix the components
    eos = preos(mix, 'qmr') # use the peng robinson equation of state

    xi = np.linspace(0, 1, 1001)

    # Make and Plot Pxy diagram
    Pmat = []
    yi = []
    for xival in xi:
        y_guess, P_guess = [xival, 1-xival], P
        yival, Pval = bubblePy(y_guess, P_guess, X=[xival, 1-xival], T=T, model=eos)
        yi.append(yival[0])
        Pmat.append(Pval)
    plt.plot(xi, Pmat)
    plt.plot(yi, Pmat)
    plt.xlabel(f'x and y of %s'%(iupacmat[i]), fontsize=16)
    plt.ylabel('pressure, P(bar)', fontsize=16)
    plt.title(f'Pxy diagram for %s + %s mixture at %s K'%(iupacmat[i],iupacmat[j],T), fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()

    # Plot Txy diagram
    Tmat = []
    yi = []
    for xival in xi:
        y_guess, T_guess = [xival, 1-xival], T
        yival, Tval = bubbleTy(y_guess, T_guess, X=[xival, 1-xival], P=P, model=eos)
        yi.append(yival[0])
        Tmat.append(Tval)
    plt.plot(xi, Tmat)
    plt.plot(yi, Tmat)
    plt.xlabel(f'x and y of %s'%(iupacmat[i]), fontsize=16)
    plt.ylabel('temperature, T(K)', fontsize=16)
    plt.title(f'Txy diagram for %s + %s mixture at %.2f bar'%(iupacmat[i],iupacmat[j],P), fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()

    # Plot xy diagram
    plt.plot(xi, yi)
    plt.plot(xi, xi)
    plt.xlabel(f'x of %s'%(iupacmat[i]), fontsize=16)
    plt.ylabel(f'y of %s'%(iupacmat[i]), fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
# %%
