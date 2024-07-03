#%% McCabe Thiele Method
from TxyPxyxy import xy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import warnings
from numpy import RankWarning, VisibleDeprecationWarning

def mccabe(comp1, comp2, xd, xb, xf = 0.5, P = None, T = None, steps = True, pointson = True):
    # Suppress RankWarning
    warnings.simplefilter('ignore', RankWarning)
    warnings.simplefilter('ignore', VisibleDeprecationWarning)
    if P is None and T is None:
        print('Please provide either a temperature or a pressure')
        return None
    if P is not None:
        T = 273.15 # K
        Pgiven = True
    elif T is not None: # elif is necessary here since T is define right before this
        P = 1e5 # Pa
        Pgiven = False
    xi, yi = xy(comp1, comp2, P = P, values = True, show = False)
    # fit a curve to the data
    z = np.polyfit(xi, yi, 30)
    p = np.poly1d(z)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(xi, p(xi))
    plt.plot(xi, xi)
    steps = 0
    x = xd
    while x > xb:
        def difference(xval):
            return p(xval) - x
        intersect = fsolve(difference, x)
        if intersect > x or intersect == x:
            print('Cannot perform McCabe-Thiele Method as equilibrium curve is below y=x at distillation composition')
            break
        # draw a horizontal line from xd,xd to the best fit curve
        plt.plot([x, intersect], [x,x], ls = '-', color = 'black') # ([initialx, finalx],[initialy, finaly])
        plt.plot([intersect, intersect], [x, intersect], ls = '-', color = 'black') # ([initialx, finalx],[initialy, finaly])
        x = intersect
        steps += 1
    #annotate at the bottom right th enumber of steps
    plt.annotate(f'steps = {steps}', xy=(0.8, 0.1), xycoords='axes fraction', fontsize=12, ha='center', va='center')
    plt.xlabel(f'Liquid mole fraction {comp1}', fontsize = 16)
    plt.ylabel(f'Vapor mole fraction {comp1}', fontsize = 16)
    plt.title(f'McCabe-Thiele Method for {comp1} + {comp2}', fontsize = 16)
    if pointson:
        #show where the xd and xb are 
        plt.plot(xd, xd, 'ro')
        plt.plot(xb, xb, 'ro')
        #label each point
        ax.annotate('xd', (xd, xd), textcoords="offset points", xytext=(0,-20), ha='center')
        ax.annotate('xb', (xb, xb), textcoords="offset points", xytext=(0,-20), ha='center')
        ax.tick_params(labelsize = 14)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
    plt.show()
    if steps:
        print(f'Steps = {steps}')

# %%
mccabe('methanol', 'water', xd = 0.95, xb = 0.1, P = 1)
# %%
