#!/usr/bin/env python
#^ This line allows running script from bash w/o python call.

#Plot omega/gamma results of linear scanscript runs. Takes in as arg final file num atm. For instance if omega_0001->omega_0015, pass in 15.

import sys
import matplotlib.pyplot as plt

kyVals = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]
omegaValsMill = [.075, .1590, .2520, .3670, .5080, .6650, .8290, .9970, 1.159, 1.314, 1.459, 1.591, 1.711, -.981, -1.054, -1.127]
gammaValsMill = [.057, .151, .285, .425, .545, .633, .686, .706, .696, .661, .607, .537, .456, 0.365, .381, .398]
omegaValsCirc = [.072,.153,.242,.352,.485,.633,.787,.942,1.091,1.234,1.487,1.365,1.597,-0.942,-1.013,-1.083]
gammaValsCirc = [.056,.146,.273,.405,.516,.596,.642,.657,.645,.610,.493,.558,.418,.356,.373,.391]

kyValsGoerler = [0.078,0.235,0.313,0.391,0.470,0.548,0.626,0.704,0.783,0.861,0.939,1.017]
omegaValsGoerler = [0.132,0.472,0.702,0.936,1.147,1.335,1.42,-0.984,-1.103,-1.221,-1.340,-1.461]
gammaValsGoerler = [0.11,0.51,0.625,0.68,0.65,0.575,0.466,0.385,0.41,0.45,0.49,0.54]

kyValsFinOld = kyVals[:-1]
kyValsFinLowV = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,.85,.9,.95,1.0,1.05,1.1,1.15,1.2]
kyValsFin = kyValsFinLowV
omegaGoerlerFinOld = [.0802,.1706,.2710,.3858,.5169,.6594,.8070,.9543,1.097,1.2318,1.3568,1.4704,1.5717,-.9493,-1.0175]
gammaGoerlerFinOld = [.0567,.1495,.277,.4106,.5244,.6070,.6555,.6722,.6615,.6281,.5767,.5117,.4371,.3390,.3516]
omegaGoerlerFinLowV = [.0802,.1705,.2707,.3853,.5164,.6589,.8068,.9544,1.0973,1.2326,1.3580,1.4720,1.5738,-.9498,-1.0181,-1.0859,-1.1531,-1.22,-1.288,-1.3583,-1.4297,-1.501,-1.5721,-1.6438]
gammaGoerlerFinLowV = [.0568,.1501,.2782,.4124,.5267,.6096,.6585,.6755,.6651,.6319,.5805,.5155,.4409,.3393,.3519,.3653,.3803,.3977,.4179,.4415,.4671,.4930,.5190,.5460]
gammaGoerlerFin = [.0553,.1463,.2736,.4066,.5184,.5987,.6449,.6598,.6478,.6136,.5619,.4967,.3396,.3563,.3732,.3909,.4101,.4309,.4534,.4783,.5057,.5353,.5670,.6009]
omegaGoerlerFin = [.0725,.1523,.2414,.3507,.4840,.6322,.7869,.9417,1.0918,1.2339,1.3665,1.4883,-.8711,-.9430,-1.0135,-1.0831,-1.1516,-1.2196,-1.2874,-1.3551,-1.4235,-1.4929,-1.5629,-1.6333]

fig = plt.figure(num=None, figsize=(12,6), dpi=100)

plt.subplot(2,1,1)
#plt.plot(kyVals, omegaValsCirc, markerfacecolor = 'green', marker = '^', color = 'green', label = 'Old Circular')
p1 = plt.plot(kyValsGoerler, omegaValsGoerler, markerfacecolor = 'black', marker = '*', color = 'black', label = 'Goerler')
p2 = plt.plot(kyValsFinLowV, omegaGoerlerFinLowV, markerfacecolor = 'red', marker = 'o', color = 'red', label = 'Low v-res')
p3 = plt.plot(kyValsFin, omegaGoerlerFin, markerfacecolor = 'yellow', marker = '.', color = 'blue', label = 'High v-res')
plt.ylabel('$\\omega$',      labelpad=14, fontsize=20)
plt.xlabel('$k_y$$\\rho_i$', labelpad=14, fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid()

plt.subplot(2,1,2)
#plt.plot(kyVals, gammaValsCirc, markerfacecolor = 'green', marker = '^', color = 'green', label = 'Old Circular')
p4 = plt.plot(kyValsGoerler, gammaValsGoerler, markerfacecolor = 'black', marker = '*', color = 'black')
p5 = plt.plot(kyValsFinLowV, gammaGoerlerFinLowV, markerfacecolor = 'red', marker = 'o', color = 'red')
p6 = plt.plot(kyValsFin, gammaGoerlerFin, markerfacecolor = 'yellow', marker = '.', color = 'blue')
plt.ylabel('$\\gamma$',      labelpad=14, fontsize=20)
plt.xlabel('$k_y$$\\rho_i$', labelpad=14, fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#Add shared legend and shift things for it.
line_labels = ['Goerler', 'Low v-res', 'High v-res']
fig.legend([p1,p2,p3],labels=line_labels,loc="upper right", fontsize=18, framealpha=1)
plt.grid()
plt.tight_layout(rect=[0,0.03,1,.95])
plt.suptitle('Cyclone ITG Mode Frequencies', fontsize=18)
plt.savefig('LinearITG_KinEl_GrowthRates.pdf')
plt.show()
