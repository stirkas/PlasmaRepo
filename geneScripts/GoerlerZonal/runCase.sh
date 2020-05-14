#!/bin/bash

#python3 ../plotPhi.py 331 192 2 5.07202 0.987922 './phi_0026_r.dat' './phi_0026_k.dat' -s -f './phiZonal3kymin.mp4'

#python3 ../plotPhi.py 248 192 2 5.07202 .592753 './phi_0027_r.dat' './phi_0027_k.dat' -s -f './phiZonal5kymin.mp4'

python3 ../plotPhi.py 82 192 2 5.07202 0.987922 './contelectrons_0009r.dat' './contelectrons_0009k.dat' -s -f './FullShearPhi.mp4'

python3 ../plotPhi.py 77 192 2 5.07202 0.987922 './contelectrons_0011r.dat' './contelectrons_0011k.dat' -s -f './StartWavePhi.mp4'

python3 ../plotPhi.py 62 192 2 5.25988 0.987922 './contelectrons_0014r.dat' './contelectrons_0014k.dat' -s -f './HalfShearPhi.mp4'

