import numpy as np
from scipy.optimize import minimize
import math

# this function holds the mass equation we want to minimize
def massFunction(x):
    areaBase = x[0]
    diam = x[1]
    lengthCyl = x[2]
    cylSpacing = x[3]
    rho = 2707
    t = 0.003

    return rho*areaBase*(t+ math.pi*lengthCyl*(diam**2)/(4.0*(cylSpacing**2)))

# this function holds the heat transfer constraint
def heatTransferFunc(x):
    areaBase = x[0]
    diam = x[1]
    lengthCyl = x[2]
    cylSpacing = x[3]
    h = 20
    k = 206
    t0 = 100
    tinf = 20

    correctLength = lengthCyl + diam/4.0

    return 50 - areaBase*((np.tanh( ((h*4.0/(k*diam))**0.5)*correctLength)/(((h*4.0/(k*diam))**0.5)*correctLength))*h*math.pi*diam*correctLength*(t0-tinf)/(cylSpacing**2)+ h*(t0-tinf)*(1-(math.pi*(diam**2))/(4.0*cylSpacing**2)))

# this function tells our program that spacing must be bigger then diameter
def constraint1(x):
    return x[3] - x[1]
# make sure that the area is bigger then spacing area
def constraint2(x):
    return x[0] - x[3]**2
# makre sure area bigger then thickness of wall
def constraint3(x):
    return x[0]-0.003

x1 = np.linspace(0.001,0.005)
y1 =
# we set up initial values
x0 = np.array([100, .01, .2, .5])
b =(0.00001, 100000.0)
bnds = (b,b,b,b)

con1 = {'type': 'eq', 'fun': heatTransferFunc}
con2 = {'type': 'ineq', 'fun': constraint1}
con3 = {'type': 'ineq', 'fun': constraint2}
con4 = {'type': 'ineq', 'fun': constraint3}
cons = [con1, con2, con3, con4]

solution_1 = minimize(massFunction, x0, method='SLSQP',bounds=bnds, constraints = cons)

print(solution_1)
