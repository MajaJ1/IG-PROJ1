def sr(data):
    '''
    Funkcja oblicza srednie współrzędne 
    '''
    Sumx = 0
    Sumy = 0
    Sumz = 0
    for sod,X,Y,Z,q in data:
        Sumx += X
        Sumy += Y
        Sumz += Z

    Sx = float(Sumx / len(data))
    Sy = float(Sumy / len(data))
    Sz = float(Sumz / len(data))

    return Sx,Sy,Sz