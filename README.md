# McCabe-Thiele Method 

Required Inputs: 

2 components, distillation composition, bottoms composition, and either a temperature or pressure, but not both

Optional Inputs: 

feed composition, murphree tray efficiency, and either reflux ratio or feed quality

Miscellaneous Inputs:

show composition point labels, show number of steps, show optimal feed location

Input Units: 

Temperature (T) in Kelvin and Pressure (P) in bar

Optional Outputs: Plot and Steps

Possible Errors: 

-If the equilibrium curve is below the y=x line at the distillation, bottoms, or the feed composition, the code will fail as it will be impossible to perform stepping off.

This will generate this error message: 

Cannot perform McCabe-Thiele Method as equilibrium curve is below y=x here
  
This may be caused by an azeotrope formation or that the first component given is less volatile than the second component
  
Consider swapping the first component with the second component as a possible solution

-If neither T nor P is inputted, the code will fail and generate this error message: Please provide either a temperature or a pressure

Dependencies:

-matplotlib.pyplot

-numpy

-scipy.optimize

-warnings

Arguments:

mccabe(comp1, comp2, xd, xb, xf = 0.5, P = None, T = None, steps = True, pointson = True)

Examples:

mccabe('methanol', 'water', xd = 0.95, xb = 0.1, P = 1)

![image](https://github.com/Chlvan/mccabepy/assets/112746859/14cdcbab-1358-4ff6-b03e-981457804715)

Steps = 4
