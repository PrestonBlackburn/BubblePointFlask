import json
import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy.optimize import fsolve


#JSONdata = '{"Temperature":74,"Hydrogen Sulfide":0,"Nitrogen":0.019264,"Carbon Dioxide":0.176679,"Methane":2.216492,"Ethane":2.131195,"Propane":2.672874,"Isobutane":1.531791,"n-Butane":3.01741,"n-Pentane":5.660723,"n-Hexane":8.777826,"n-Heptane":73.795745,"n-Octane":0,"n-Nonane":0,"n-Decane":0,"Benzene":0,"Toluene":0,"Xylenes":0,"Cyclohexane":0}'
JSONdata = '{"Benzene": 0,"Carbon Dioxide": 0.176679,"Cyclohexane": 0,"Ethane": 2.131195,"Hydrogen Sulfide": 0,"Isobutane": 1.531791,"Methane": 2.216492,"Nitrogen": 0.019264,"Propane": 2.672874,"Temperature": 74,"Toluene": 0,"Xylenes": 0,"n-Butane": 3.01741,"n-Decane": 0,"n-Heptane": 73.795745,"n-Hexane": 8.777826,"n-Nonane": 0,"n-Octane": 0,"n-Pentane": 5.660723}'


def main(JSONdata):
    (dfIn,Tfeed) = dataFormat(JSONdata)
    bubble = Bubble(dfIn,Tfeed)
    bubblePoint = bubble[0]*14.5038-14.73
    print(dfIn)
    print(Tfeed)
    print(bubblePoint)

def dataFormat(JSONdata):
    data = json.loads(JSONdata)

    T_feed = data['Temperature']
    Tfeed = (T_feed-32)*5/9+273.15

    composition = {'Component':['Hydrogen Sulfide', 'Nitrogen', 'Carbon Dioxide', 'Methane', 'Ethane', 'Propane', 'Isobutane', 'n-Butane', 'n-Pentane', 'n-Hexane', 'n-Heptane', 'n-Octane', 'n-Nonane', 'n-Decane', 'Benzene', 'Toluene', 'Xylenes', 'Cyclohexane']}

    dummy = []
    for component in composition['Component']:
        dummy.append(data[component])

    molvals = {'Mol%':dummy}
    composition.update(molvals)  

    dfIn = pd.DataFrame(composition)
    dfIn['MolFrac']= dfIn['Mol%']/(100)

    return dfIn, Tfeed

def Bubble(dfIn,T):
    constants = Constants()
    P1 = dfIn['MolFrac']*np.exp(np.log(constants.df['criticalPress']) + np.log(10)*(7/3)*(1 + constants.df['accentricFac'])*(1-constants.df['criticalTemp']/T))
    P1 = P1.sum()
    guess = P1
 
    data = (dfIn, T)
    ro = fsolve(BubbleFunctionToOptimize,guess, args = data)

    return(ro)


def BubbleFunctionToOptimize(BubbleGuess,*data):
    (dfIn, T) = data
    constants = Constants()
    initialize = bubble_initialize_function(dfIn,constants,BubbleGuess,T)
    
    mixture = mixture_function(dfIn['MolFrac'],initialize['df'],constants,BubbleGuess,T)
    Z = cubic_solver_function(mixture['A'], mixture['B'])
    ZL = min(Z)
    fugacity_init = fugacity_initialize_function(initialize['df'],constants,BubbleGuess,T)
    fugacityA = np.dot(fugacity_init['fugacityMatrix'],dfIn['MolFrac'])
    phiL = fugacity_function(fugacityA,fugacity_init['fugacityb'],ZL,mixture['A'], mixture['B'])

    Yraw_old = initialize['sumY']
    yold = initialize['df']['Y']
    max_iter = 50
    i = 0
    found = False
    while  not found and i < max_iter:
        mixture = mixture_function(yold,initialize['df'],constants,BubbleGuess,T)
        Z = cubic_solver_function(mixture['A'], mixture['B'])
        ZV = max(Z)
        fugacityA = np.dot(fugacity_init['fugacityMatrix'],yold)
        phiV = fugacity_function(fugacityA,fugacity_init['fugacityb'],ZV, mixture['A'], mixture['B'])
        Yraw = dfIn['MolFrac']*phiL/phiV
        Yraw_sum = Yraw.sum()
        diff = abs(Yraw_sum-Yraw_old)
        if diff <= 1.0e-16:
            found = True
        Yraw_old = Yraw_sum
        Y = Yraw/Yraw_sum
        yold = Y
        i = i+1
    if found:
        return (Yraw_sum-1)
    else:
        return np.nan

class Constants:
    def __init__(self):
        constants = {'Component':['Hydrogen Sulfide', 'Nitrogen', 'Carbon Dioxide', 'Methane', 'Ethane', 'Propane', 'Isobutane', 'n-Butane', 'n-Pentane', 'n-Hexane', 'n-Heptane', 'n-Octane', 'n-Nonane', 'n-Decane', 'Benzene', 'Toluene', 'Xylenes', 'Cyclohexane'],
                    'antoineA':[16.104, 14.9342, 22.5898, 15.2243, 15.6637, 15.726, 15.5381, 15.6782, 15.8333, 15.8366, 15.8737, 15.9426, 15.9671, 16.0114, 15.9008, 16.0137, 16.1156, 15.7527],
                    'antoineB':[1768.69, 588.72, 3103.39, 897.84, 1511.42, 1872.46, 2032.73, 2154.9, 2477.07, 2697.55, 2911.32, 3120.29, 3291.45, 3456.8, 2788.51, 3096.52, 3395.57, 2766.63],
                    'antoineC':[-26.06, -6.6, -0.16, -7.16, -17.16, -25.16, -33.15, -34.42, -39.94, -48.78, -56.51, -63.63, -71.33, -78.67, -52.36, -53.67, -59.46, -50.5],
                    'molecularWeight':[34.0809, 28.0134, 44.0095, 16.0425, 30.069, 44.0956, 58.1222, 58.1222, 72.1488, 86.1754, 100.2019, 114.2285, 128.2551, 142.2817, 78.1118, 92.1384, 106.165, 84.1595],
                    'criticalPress':[89.36865, 33.943875, 73.7646, 46.00155, 48.83865, 42.455175, 36.477, 37.996875, 33.741225, 29.688225, 27.35775, 24.824625, 23.1021, 21.0756, 48.939975, 41.13795, 37.2876, 40.73265],
                    'criticalTemp':[373.2, 126.2, 304.2, 190.6, 305.4, 369.8, 408.1, 425.2, 469.6, 507.4, 540.2, 568.8, 594.6, 617.6, 562.1, 591.7, 630.2, 553.4],
                    'accentricFac':[0.1, 0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.251, 0.296, 0.351, 0.394, 0.444, 0.49, 0.212, 0.257, 0.314, 0.213],
                    'criticalCompress':[0.284, 0.29, 0.274, 0.288, 0.285, 0.281, 0.283, 0.274, 0.262, 0.26, 0.263, 0.259, 0.26, 0.247, 0.271, 0.264, 0.263, 0.273]}
        self.df = pd.DataFrame(constants)

        self.interact = np.array([[0, 0.1767, 0.0974, 0.084, 0.0833, 0.0878, 0.0474, 0.0006, 0.063, 0.0032, 0.0046, 0.0059, 0.0072, 0.0333, 0.0057, 0.0072, 0.0092, 0.0053],
                            [0.1767, 0, -0.017, 0.0311, 0.0515, 0.0852, 0.1033, 0.0711, 0.01, 0.0603, 0.1441, -0.41, 0.0722, 0.1122, 0.1641, 0.0726, 0.0778, 0.0682],
                            [0.0974, -0.017, 0, 0.0919, 0.1322, 0.1241, 0.12, 0.1333, 0.0063, 0.0087, 0.1, 0.0127, 0.0145, 0.1141, 0.0774, 0.1056, 0.0173, 0.1052],
                            [0.084, 0.0311, 0.0919, 0, -0.0026, 0.014, 0.0256, 0.0133, 0.023, 0.0313, 0.0352, 0.0496, 0.0474, 0.0411, 0.0363, 0.097, 0.0454, 0.0389],
                            [0.0833, 0.0515, 0.1322, -0.0026, 0, 0.0013, 0.0029, 0.0038, 0.0063, 0.0086, 0.0067, 0.0185, 0.0145, 0.0144, 0.0322, 0.0145, 0.0173, 0.0178],
                            [0.0878, 0.0852, 0.1241, 0.014, 0.0013, 0, 0.0003, 0.0033, 0.0019, 0.0034, 0.0056, 0.0061, 0.0074, 0.0086, 0.0233, 0.0074, 0.0094, 0.0055],
                            [0.0474, 0.1033, 0.12, 0.0256, 0.0029, 0.0003, 0, -0.0004, 0.0007, 0.0016, 0.0026, 0.0037, 0.0047, 0.0057, 0.0035, 0.0047, 0.0063, 0.0032],
                            [0.0006, 0.0711, 0.1333, 0.0133, 0.0038, 0.0033, -0.0004, 0, 0.0003, 0.0011, 0.0033, 0.0074, 0.0037, 0.0078, 0.0026, 0.0037, 0.0052, 0.0024],
                            [0.063, 0.01, 0.0063, 0.023, 0.0063, 0.0019, 0.0007, 0.0003, 0, 0.0002, 0.0074, 0.0012, 0.0018, 0.0024, 0.0174, 0.0018, 0.0029, 0.0037],
                            [0.0032, 0.0603, 0.0087, 0.0313, 0.0086, 0.0034, 0.0016, 0.0011, 0.0002, 0, -0.0078, 0.0004, 0.0008, 0.0012, 0.0093, 0.0008, 0.0015, -0.003],
                            [0.0046, 0.1441, 0.1, 0.0352, 0.0067, 0.0056, 0.0026, 0.0033, 0.0074, -0.0078, 0, 0.0001, 0.0003, 0.0006, 0.0011, 0.0003, 0.0008, 0],
                            [0.0059, -0.41, 0.0127, 0.0496, 0.0185, 0.0061, 0.0037, 0.0074, 0.0012, 0.0004, 0.0001, 0, 0.0001, 0.0002, 0.003, 0.0001, 0.0003, 0],
                            [0.0072, 0.0722, 0.0145, 0.0474, 0.0145, 0.0074, 0.0047, 0.0037, 0.0018, 0.0008, 0.0003, 0.0001, 0, 0, 0.0001, 0, 0.0001, 0.0002],
                            [0.0333, 0.1122, 0.1141, 0.0411, 0.0144, 0.0086, 0.0057, 0.0078, 0.0024, 0.0012, 0.0006, 0.0002, 0, 0, 0.0003, 0.0001, 0, 0.0004],
                            [0.0057, 0.1641, 0.0774, 0.0363, 0.0322, 0.0233, 0.0035, 0.0026, 0.0174, 0.0093, 0.0011, 0.003, 0.0001, 0.0003, 0, 0.0001, 0.0004, 0.0126],
                            [0.0072, 0.0726, 0.1056, 0.097, 0.0145, 0.0074, 0.0047, 0.0037, 0.0018, 0.0008, 0.0003, 0.0001, 0, 0.0001, 0.0001, 0, 0.0001, 0.0001],
                            [0.0092, 0.0778, 0.0173, 0.0454, 0.0173, 0.0094, 0.0063, 0.0052, 0.0029, 0.0015, 0.0008, 0.0003, 0.0001, 0, 0.0004, 0.0001, 0, 0.0006],
                            [0.0053, 0.0682, 0.1052, 0.0389, 0.0178, 0.0055, 0.0032, 0.0024, 0.0037, -0.003, 0, 0, 0.0002, 0.0004, 0.0126, 0.0001, 0.0006, 0]])
        self.R = 8.3145


def bubble_initialize_function(dfIn,constants,BubbleGuess,T):
    initialize = {'Component':['Hydrogen Sulfide', 'Nitrogen', 'Carbon Dioxide', 'Methane', 'Ethane', 'Propane', 'Isobutane', 'n-Butane', 'n-Pentane', 'n-Hexane', 'n-Heptane', 'n-Octane', 'n-Nonane', 'n-Decane', 'Benzene', 'Toluene', 'Xylenes', 'Cyclohexane']}
    dfInit = pd.DataFrame(initialize)
    dfInit['P1'] = dfIn['MolFrac']*np.exp(np.log(constants.df['criticalPress']) + np.log(10)*(7/3)*(1 + constants.df['accentricFac'])*(1-constants.df['criticalTemp']/T))
    dfInit['K'] = np.exp(np.log(constants.df['criticalPress']/BubbleGuess) + np.log(10)*(7/3)*(1 + constants.df['accentricFac'])*(1-constants.df['criticalTemp']/T))
    dfInit['Yint'] = dfInit['K'] * dfIn['MolFrac']
    sumY = dfInit['Yint'].sum(axis = 0)
    dfInit['Y'] = dfInit['Yint']/sumY
    dfInit['κ'] = 0.37464 + 1.54226*constants.df['accentricFac'] - 0.26992*constants.df['accentricFac'].pow(2)
    dfInit['α'] = (1 + dfInit['κ'] * (1 - (T/constants.df['criticalTemp']).pow(0.5))).pow(2)
    dfInit['a'] = 0.45724 * dfInit['α'] * (constants.R * constants.df['criticalTemp']).pow(2)/(constants.df['criticalPress'] * 100000)
    dfInit['b'] = 0.07780 * constants.R * constants.df['criticalTemp'] / (constants.df['criticalPress'] * 100000)

    d = dict()
    d['df'] = dfInit
    d['sumY'] = sumY
    return d


def mixture_function(molfrac,dfInit,constants,BubbleGuess,T):
    # aij = [(ai.aj)0.5(1 - kij)] = aji
    mixture = np.outer(dfInit['a'],dfInit['a'])
    mixture = np.power(mixture,0.5)
    mixture = np.multiply(mixture, np.subtract(1,constants.interact))
    mixture = np.multiply(mixture,np.outer(molfrac,molfrac))
    asum  = mixture.sum()
    b = dfInit['b']*molfrac
    sumb = b.sum()

    A = asum * BubbleGuess * 100000 / (constants.R * T) ** 2
    B = sumb * BubbleGuess * 100000 / (constants.R * T)

    d = dict()
    d['A'] = A
    d['B'] = B
    return d


def cubic_solver_function(A,B):
    C2 = B - 1
    C1 = A-3*B**2-2*B
    C0 = -1*A*B+B**2+B**3

    Q1 = C2*C1/6 - C0/2 - C2 ** 3/27
    P1 = C2**2/9-C1/3
    D = Q1**2-P1**3 

    if D >= 0:
        Z = [(Q1+D**0.5)**(1/3) + pow(abs(Q1-D**0.5),float(1)/3) *np.sign(Q1-D**0.5) - C2/3, np.nan, np.nan]
    else:
        t1 = Q1**2/P1**3
        t2 = (1-t1)**0.5/t1**0.5*Q1/abs(Q1)
        phi = np.arctan(t2)
        if phi<0:
            phi = phi + np.pi
        Z0 = 2*P1**0.5*np.cos(phi/3)-C2/3
        Z1 = 2*P1**0.5*np.cos((phi+2*np.pi)/3)-C2/3
        Z2 = 2*P1**0.5*np.cos((phi+4*np.pi)/3)-C2/3
        Z = [Z0,Z1,Z2]
    
    return(Z)


def fugacity_initialize_function(dfInit,constants,BubbleGuess,T):
    fugacity = np.outer(dfInit['a'],dfInit['a'])
    fugacity = np.power(fugacity,0.5)
    fugacity = np.multiply(fugacity, np.subtract(1,constants.interact))
    fugacity = fugacity*BubbleGuess*100000/(constants.R*T)**2
    fugacityb = dfInit['b']*BubbleGuess*100000/(constants.R*T)

    d = dict()
    d['fugacityMatrix'] = fugacity
    d['fugacityb'] = fugacityb
    return d


def fugacity_function(fugacityA,fugacityb,Z,A,B):
    V1  = fugacityb/B
    V2 = Z - 1
    V3 = np.log(Z - B)
    V4 = A*(2*fugacityA/A-fugacityb/B)
    V5 = Z + 2.41421536*B
    V6 = Z - 0.41421536*B
    V7 = 2.82842713*B
    phiL = np.exp(V1*V2-V3-V4*np.log(V5/V6)/V7)
    return phiL

main(JSONdata)