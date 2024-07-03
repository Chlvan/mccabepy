#%% McCabe Thiele Method
from TxyPxyxy import xy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import warnings
from numpy import RankWarning, VisibleDeprecationWarning

def mccabe(comp1, comp2, xd, xb, xf = 0.5, P = None, T = None, R = None, q = None, stages = True, feedstage = True, pointson = True, showplot = True, values = False):

    # Suppress Errors
    warnings.simplefilter('ignore', RankWarning)
    warnings.simplefilter('ignore', VisibleDeprecationWarning)
    warnings.simplefilter('ignore', RuntimeWarning)

    # Check if P or T is given
    if P is None and T is None:
        print('Please provide either a temperature or a pressure')
        return None
    if P is not None:
        T = 273.15 # K
        Pgiven = True
    elif T is not None: # elif is necessary here since T is define right before this
        P = 1e5 # Pa
        Pgiven = False

    #Check for q and R value condition
    if q is not None and R is None:
        print('Please provide a R value')
        return None
    if R is not None and q is None:
        print('Please provide a q value')
        return None
    if q is not None and R is not None:
        if q == 1: # this will make the slopes infinity
            feedslope = 1e10
        else:
            feedslope = q/(q-1)
        if R == -1: # this will make the slopes infinity
            rectifyslope = 1e10
        else:
            rectifyslope = R/(R + 1)
        if feedslope == rectifyslope:
            print('Please provide valid q and R values')
            return None

    # Get the equilibrium data
    xi, yi = xy(comp1, comp2, P = P, values = True, show = False)

    # Fit a curve to the data
    z = np.polyfit(xi, yi, 30)
    p = np.poly1d(z)

    # Initialize the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(xi, p(xi))
    plt.plot(xi, xi)

    # Define and plot the feed, rectifying, and stripping lines and intersection
    if q is not None and R is not None:
        def rectifying(xval):
            return R/(R+1)*xval + xd/(R+1)
        def feed(xval):
            return q/(q-1)*xval - xf/(q-1)
        def feedrectintersection(xval):
            return q/(q-1)*xval - xf/(q-1) - R/(R+1)*xval - xd/(R+1)
        xguess = xf*q # q>1 means xsol > xf, and q<1 means xsol < xf
        if R == -1:
            xsol = xd
            ysol = feed(xsol)
        elif q == 1:
            xsol = xf
            ysol = rectifying(xsol)
        else:
            xsol = fsolve(feedrectintersection, xguess)  # fsolve returns an array, get the first element
            ysol = rectifying(xsol)
        # define the stripping line as the line from xb to the intersection of the rectifying line and the feed line
        def stripping(xval):
            return (ysol-xb)*(xval-xb)/(xsol-xb)+xb
        plt.plot([xf, xsol], [xf, ysol], ls=':', color='black', label='Feed to Rectifying Line')
        plt.plot([xd, xsol], [xd, ysol], ls='--', color='black', label='Distillate to Feed Line')
        plt.plot([xb, xsol], [xb, ysol], ls='-.', color='black', label='Bottom to Feed Line')

    # Perform Step-off process for no feed condition
    stages = 0 # initialize the number of stages
    x = xd # start off at the distillation composition
    if q is None and R is None:
        xs = []
        ys = []
        xs.append(x)
        ys.append(x)
        while x > xb:
            def difference(xval):
                return p(xval) - x
            intersect = fsolve(difference, 0)
            if intersect > x or intersect == x:
                print('Cannot perform McCabe-Thiele Method as equilibrium curve is below y=x at distillation composition')
                break
            # draw a horizontal line from y=x to the best fit curve
            plt.plot([x, intersect], [x,x], ls = '-', color = 'black') # ([initialx, finalx],[initialy, finaly])
            plt.plot([intersect, intersect], [x, intersect], ls = '-', color = 'black') # ([initialx, finalx],[initialy, finaly])
            x = intersect
            xs.append(x)
            ys.append(x)
            stages += 1

    # Perform step off for feed condition
    else:
        y = xd
        xs = []
        ys = []
        xs.append(x)
        ys.append(y)
        feedstage = 1
        while x > xb:
            def difference(xval):
                return p(xval) - y
            intersect = fsolve(difference, 0)
            if intersect > x or intersect == x:
                print('Cannot perform McCabe-Thiele Method as equilibrium curve is below y=x at distillation composition')
                break
            plt.plot([x, intersect], [y,y], ls = '-', color = 'black')
            if intersect > xsol:
                # draw a horizontal line from rectifying line to the best fit curve 
                plt.plot([intersect, intersect], [y, rectifying(intersect)], ls = '-', color = 'black')
                x = intersect
                y = rectifying(intersect)
                xs.append(x)
                ys.append(x)
                stages += 1
                feedstage += 1
            else:
                # draw a horizontal line from stripping line to the best fit curve 
                plt.plot([intersect, intersect], [y, stripping(intersect)], ls = '-', color = 'black')
                x = intersect
                y = stripping(intersect)
                xs.append(x)
                ys.append(x)
                stages += 1

    # Add labels to the plot
    plt.xlabel(f'Liquid mole fraction {comp1}', fontsize = 16)
    plt.ylabel(f'Vapor mole fraction {comp1}', fontsize = 16)
    plt.title(f'McCabe-Thiele Method for {comp1} + {comp2}', fontsize = 16)

    # Allow for points to be shown
    if pointson:
        plt.plot(xd, xd, 'ro')
        plt.plot(xb, xb, 'ro')
        if q is not None and R is not None:
            plt.plot(xf, xf, 'ro')
        ax.annotate('xd', (xd, xd), textcoords="offset points", xytext=(0,-15), ha='center')
        if xb < 0.05:
            ax.annotate('xb', (xb, xb), textcoords="offset points", xytext=(15,-3), ha='center')
        else:
            ax.annotate('xb', (xb, xb), textcoords="offset points", xytext=(0,-15), ha='center')
        if q is not None and R is not None:
            ax.annotate('xf', (xf, xf), textcoords="offset points", xytext=(0,-15), ha='center')
        ax.tick_params(labelsize = 14)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)

    # Creating a nice output boxes
    textstr1 = '\n'.join((
        r'$Input:$',
        r'$x_D=%.2f$' % (xd, ),
        r'$x_F=%.2f$' % (xf, ),
        r'$x_B=%.2f$' % (xb, ),
        ))  # Close the parentheses here
    if q is not None and R is not None:
        textstr1 += '\n' + '\n'.join((
            r'$q=%.2f$' % (q, ),
            r'$R=%.2f$' % (R, ),
        ))  # Close the parentheses here
    textstr2 = '\n'.join((
        r'$Output:$',
        r'$Stages=%.0f$' % (stages, ),
        ))
    if q is not None and R is not None:
        textstr2 += '\n' + '\n'.join((
            r'$Feed \:Stage=%.0f$' % (feedstage, ),
        ))
    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='green', alpha=0.2)
    plt.text(0.55, 0.05, textstr1, fontsize=10, verticalalignment='bottom', bbox=props)
    plt.text(0.75, 0.05, textstr2, fontsize=10, verticalalignment='bottom', bbox=props)

    # Show the plot
    if showplot:
        plt.show()
    
    # Return the values
    if values:
        if q is not None and R is not None:
            return xs, ys, stages, feedstage, xf, xsol, ysol
        else:
            return xs, ys

# %%
mccabe('methanol', 'water', xd = 0.95, xb = 0.01, P = 1, q = 1, R = 3)
mccabe('methanol', 'water', xd = 0.95, xb = 0.1, P = 1)
# %%
