import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

"""Script to produce model results for a stirred mixed reactor"""

plt.ion()
sns.set_context('poster',font_scale=1.25)
plt.rcParams['font.family'] = 'Cambria'
plt.rcParams['svg.fonttype'] = 'none'

tRes = np.logspace(-6,6,130)
# Stirred reactor at steady state
k = 2000  # 1/s/M
cB0 = 3E-3  # M
cA0 = 0.5E-3  # M
cA = cA0/(1+k*cB0*tRes)  # Normalized to
Da = k*cA0*tRes
r = k*cB0*cA

f1, ax1 = plt.subplots(1,1,figsize=(12,10))
#ax1.plot(tRes, cA/cA0, label='Reactant Concentration')
ax1.plot(Da,k*cB0*cA/(k*cB0*cA0), label='Reaction Rate')
ax1.plot(Da, 1-cA/cA0, label = 'Product')
plt.xlabel('Dahmkohler Number')
plt.ylabel('Normalized Value (-)')
plt.xscale('log')
plt.legend()
sns.despine(f1)

plt.show()