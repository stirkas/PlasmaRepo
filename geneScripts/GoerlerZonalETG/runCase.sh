#!/bin/bash

#Ad el real/k space.
#python3 ../plotPhi.py 637 192 16 5.63558 2.96377 './adElPhi.dat' './adElPhik.dat' -s -f './phiZonalETG.mp4
#Ad el k-space + flux.
python3 ../plotPhi.py 1879 192 16 5.63558 2.96377 './adIonPhi.dat' './adIonPhik.dat' -q ./adIonFlux.dat -s -f './phiZonalETGwFluxAdIon.mp4'

#Kin el real/k space.
#Kin el k-space + flux.
python3 ../plotPhi.py 1879 192 16 5.63558 2.96377 './kinIonPhi.dat' './kinIonPhik.dat' -q ./kinIonFlux.dat -s -f './phiZonalETGwFluxKinIon.mp4'
