import math
import numpy as np

def l_NEU(fi, l, X_sr, Y_sr, Z_sr, data):
    '''
    funkcja przelicza z współrzedne 
    do układu  (East (E), North (N), Up (U))
    '''
    X=[]
    Y=[]
    Z=[]
    for line in data:
        X.append(line[1])
        Y.append(line[2])
        Z.append(line[3])
    X_sr = float(X_sr)
    Y_sr = float(Y_sr)
    Z_sr = float(Z_sr)
    dX = []
    dY = []
    dZ = []
    NEU = []
    dN = []
    dE = []
    dU = []
    for x, y, z in zip(X, Y, Z):
        delta_X = x - X_sr
        delta_Y = y - Y_sr
        delta_Z = z - Z_sr
        dX.append(delta_X)
        dY.append(delta_Y)
        dZ.append(delta_Z)
    Rt = np.matrix([((-math.sin(fi) * math.cos(l)), (-math.sin(fi) * math.sin(l)), (math.cos(fi))),
                    ((-math.sin(l)), (math.cos(l)), (0)),
                    ((math.cos(fi) * math.cos(l)), (math.cos(fi) * math.sin(l)), (math.sin(fi)))])
    for x, y, z in zip(dX, dY, dZ):
        d = np.matrix([x, y, z])
        d = d.T
        neu = Rt * d
        NEU.append(neu)
        dN.append(float(neu[0]))
        dE.append(float(neu[1]))
        dU.append(float(neu[2]))

    return (NEU, dN, dE, dU)