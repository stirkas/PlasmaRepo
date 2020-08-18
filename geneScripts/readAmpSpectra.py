#!/usr/bin/env python

def readAmpSpectra(fileName):
   f = open(fileName, 'r')

   krho  = []
   phi   = []
   n     = []
   Tpar  = []
   Tperp = []
   qpar  = []
   qperp = []
   upar  = []
   allSpectraData = []

   for line in f.readlines():
      if (line[0] != '#'): #Ignore comment lines.
         if (len(line.split()) > 0): #Ignore empty lines.
            allSpectraData.append(line.split())
   f.close()
    
   for data in allSpectraData:
      krho.append(float(data[0]))
      phi.append(float(data[1]))
      n.append(float(data[2]))
      Tpar.append(float(data[3]))
      Tperp.append(float(data[4]))
      qpar.append(float(data[5]))
      qperp.append(float(data[6]))
      upar.append(float(data[7]))
   
   #Find index splitting ky data and kx data.
   index = 0
   for i in range(len(krho)):
      if (krho[i+1] < krho[i]):
         index = i+1
         break

   kySpectra = [krho[:index], phi[:index], n[:index], Tpar[:index], Tperp[:index], qpar[:index], qperp[:index], upar[:index]]
   kxSpectra = [krho[index:], phi[index:], n[index:], Tpar[index:], Tperp[index:], qpar[index:], qperp[index:], upar[index:]]

   return [kxSpectra, kySpectra]