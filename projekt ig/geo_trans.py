# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:19:15 2022

@author: Maja
"""
import numpy as np
import math
from math import sin, cos, sqrt, atan, atan2, degrees, radians, tan
dane = np.array([[3664940.500,1409153.590,5009571.170],
               [3664940.510,1409153.580,5009571.167],
               [3664940.520,1409153.570,5009571.167],
               [3664940.530,1409153.560,5009571.168],
               [3664940.520,1409153.590,5009571.170],
               [3664940.514,1409153.584,5009571.166],
               [3664940.525,1409153.575,5009571.166],
               [3664940.533,1409153.564,5009571.169],
               [3664940.515,1409153.590,5009571.170],
               [3664940.514,1409153.584,5009571.169],
               [3664940.515,1409153.595,5009571.169],
               [3664940.513,1409153.584,5009571.171]])


plik = "wsp_inp.txt"

tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)


class Transformacje():
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        global a, b, ecc, ecc2

    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            

    def wsp_geod2XYZ(self, fi, lam, h, a, ecc2):
        N = self.a/math.sqrt(1-self.ecc2 * math.sin(fi)**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N*(1-self.ecc2) + h) * math.sin(fi)
        return (X, Y, Z)
        
        #X, Y, Z = wsp_geod2XYZ(fi, lam, h, a, self.ecc2)
        
        print('X =', round(X, 3), 'Y =', round(Y, 3), 'Z =', round(Z, 3))
            
    def uklad_1992(self, fi, lam, m_0, output='dec_degree'):
        N = self.a/(sqrt(1-self.ecc2 * sin(fi)**2))
        t = tan(fi)
        n2 = self.ecc2/(1-self.ecc2)* cos(lam)**2
        lam_0 = radians(19)
        l = lam - lam_0
        
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
        
        sigma = self.a* ((A_0*fi) - (A_2*sin(2*fi)) + (A_4*sin(4*fi)) - (A_6*sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        return x92, y92 
        
    
    
    def u2(self,fi, lam, m_0, output='dec_degree'): 
   
        N = self.a/(math.sqrt(1-self.ecc2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.ecc2 * np.cos(lam)**2
        lam = math.degrees(lam)
        
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24  
            
        #lam = math.radians(lam)
        #lam_0 = math.radians(lam_0)
        #l = lam - lam_0
        
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
        
        l =  math.radians(lam) - math.radians(lam_0)

        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
        
        return x00, y00 


#x, y = uklad_2000(fi_A, lam_A, a, self.ecc2, m_0)
#x1, y1 = uklad_2000(fi_B, lam_B, a, self.ecc2, m_0) 


a = Transformacje()
print(type(float(dane[1,0])))
i = 0
wsp_hirv = []

for i in dane:
    h = a.xyz2plh(i[0], i[1], i[2])
    wsp_hirv.append(h)
    print('wsp_hirv = ', wsp_hirv)
        
    
a = Transformacje()
print(type(float(dane[1,0])))
i = 0
u_1992 = []

for i in dane:
    j = a.uklad_1992(i[0], i[1], i[2])
    u_1992.append(j)
    print('u_1992 = ', u_1992)
    

#a = Transformacje()
#print(type(float(dane[1,0])))
#i = 0
#u2 = []

#for i in dane:
 #   j = a.u2(i[0], i[1], i[2])
  #  u2.append(j)
   # print('u2 = ', u2)
        

np.savetxt("wsp_out.txt", tablica, delimiter=',', fmt = ['%10.2f', '%10.2f', '%10.3f'], header = 'konwersja współrzednych geodezyjnych \\ ')




# if __name__ == "__main__":
#     # utworzenie obiektu
#     geo = Transformacje(model = "wgs84")
#     # dane XYZ geocentryczne
#     X = 3664940.500; Y = 1409153.590; Z = 5009571.170
#     phi, lam, h = geo.xyz2plh(X, Y, Z)
#     print(phi, lam, h)
#     #phi, lam, h = geo.xyz2plh2(X, Y, Z)
#     #print(phi, lam, h)
