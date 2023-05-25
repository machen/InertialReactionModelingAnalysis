import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

"""Script to produce model results for a stirred mixed reactor"""

plt.ion()

tRes = np.logspace(-6,6,130)
# Stirred reactor at steady state
k = 2000  # 1/s/M
cB0 = 3E-3  # M
cA0 = 0.5E-3  # M
cA = cA0/(1+k*cB0*tRes)  # Normalized to
Da = k*cA0*tRes
r = k*cB0*cA

f1, ax1 = plt.subplots(1,1)
ax1.plot(tRes, cA/cA0)
plt.xlabel('Residence Time (sec)')
plt.ylabel('Normalized Reactant concentration (-)')
plt.xscale('log')
sns.despine(f1)

f2,ax2 = plt.subplots(1,1)
ax2.plot(Da,k*cB0*cA/(k*cB0*cA0))
plt.xlabel('Dahmkohler Number')
plt.ylabel('Normalized Reaction Rate (-)')
plt.xscale('log')
sns.despine(f2)

plt.show()