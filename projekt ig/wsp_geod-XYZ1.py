import math

#Zmiana fi lam h na XYZ

a= 6378137
e2= 0.00669438002290
fi = 0.9219852924569152
lam = 0.25014684307590096
h = 200

def wsp_geod2XYZ(fi, lam, h, a, e2):
    N = a/math.sqrt(1-e2*math.sin(fi)**2)
    X = (N + h) * math.cos(fi) * math.cos(lam)
    Y = (N + h) * math.cos(fi) * math.sin(lam)
    Z = (N*(1-e2) + h) * math.sin(fi)
    return(X, Y, Z)

X, Y, Z = wsp_geod2XYZ(fi, lam, h, a, e2)

print('X =', round(X, 3), 'Y =', round(Y, 3), 'Z =', round(Z, 3))