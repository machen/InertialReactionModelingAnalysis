import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

dataFile = "..\\SimulationSummaryData_v6.xlsx"
data = pd.read_excel(dataFile, usecols="A:N,T:Y",
                     names=["FileName", "r1", "r2", "d", "k",
                            "c", "Re", "ReP", "FlowCond", "uInlet", "EstRT",
                            "MRT", "MRTerr", "TimeRatio",
                            "scalDisTCPO", "scalDisH2O2",
                            "scalDisProd", "scalDisConserv",
                            "dCdtMean", "dCdtStd"],
                     skiprows=2, engine="openpyxl")
data.sort_values(by=["FlowCond","ReP"], inplace=True)

plt.ion()
sns.set_context('paper',font_scale=1.25)
plt.rcParams['font.family'] = 'Cambria'
plt.rcParams['svg.fonttype'] = 'none'

# Figure 2: Mean dCdt results
f2, (ax2a, ax2b) = plt.subplots(ncols=1, nrows=2, figsize=(5,10))
subData100 = data.loc[(data.d == 100) & (data.k==2000) & (data.FlowCond=='NS'), :]
ax2a.plot(subData100.ReP,
          subData100.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, label="Sim 100um")
ax2b.plot(subData100.ReP,
          subData100.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, label="Sim 100um")
subData25 = data.loc[(data.d == 25) & (data.k==2000) & (data.FlowCond=='NS'), :]
ax2b.plot(subData25.ReP,
          subData25.dCdtMean/max(subData25.dCdtMean), ls='-',
          marker=None, label="Sim 25um")
subData100Stokes = data.loc[(data.d == 100) & (data.k==2000) & (data.FlowCond=='Stokes'), :]
ax2a.plot(subData100Stokes.ReP,
          subData100Stokes.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, label="Sim 100um")
ax2b.plot(subData100Stokes.ReP,
          subData100Stokes.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, label="Sim 100um")
# subData25Stokes = data.loc[(data.d == 25) & (data.k==2000) & (data.FlowCond=='Stokes'), :]
# ax2b.plot(subData25.ReP,
#           subData25.dCdtMean/max(subData25.dCdtMean), ls='-',
#           marker=None, label="Sim 25um")


ax2a.legend()
ax2a.set_xlabel('Reynolds Number')
ax2a.set_ylabel('Mean Reaction Rate Normalized to Max (.)')
ax2b.legend()
ax2b.set_xlabel('Reynolds Number')
ax2b.set_ylabel('Mean Reaction Rate Normalized to Max (.)')
ax2a.set_xlim([-1, 105])
ax2a.set_ylim([0.3, 1.05])
ax2b.set_xlim([-1, 105])
ax2b.set_ylim([0.3, 1.05])

# Figure 3: MRT and Scalar Dissipation
f3, (ax3a, ax3b) = plt.subplots(ncols=1, nrows=2, figsize=(5,10))
tReact = 1/3E-3/2000  # Half life assuming well mixed and equal reactants
ax3c = ax3b.twinx()
ax3a.plot(subData100.ReP, subData100.MRT, marker="o",
          ls='none', label="100 um Gap")
ax3a.plot(subData25.ReP, subData25.MRT, marker="o",
          ls='none', label="25 um Gap")
ax3a.plot(subData100Stokes.ReP, subData100Stokes.MRT, marker="o",
          ls='none', label="100 um Gap Stokes")
# ax3a.plot(subData25Stokes.ReP, subData25Stokes.MRT, marker="o",
#           ls='none', label="25 um Gap")
ax3a.plot([subData100.ReP.min(), subData100.ReP.max()],[tReact, tReact],
          ls='--', label="Half life Equal Mixed")
ax3a.plot([subData100.ReP.min(), subData100.ReP.max()],[9*tReact, 9*tReact],
          ls='--', label="90% Completion Equal Mixed")
ax3a.set_xlabel('Reynolds Number')
ax3a.set_ylabel('Mean residence time (s)')
ax3a.set_yscale('log')
ax3a.set_ylim([0.0005, 100])
ax3a.legend()

ax3b.plot(subData100.ReP, subData100.scalDisConserv,
          marker="o", ls='none', label="100 um Gap")
ax3b.plot(subData25.ReP, subData25.scalDisConserv,
          marker="o", ls='none', label="25 um Gap")
ax3b.plot(subData100Stokes.ReP, subData100Stokes.scalDisConserv,
          marker="o", ls='none', label="100 um Gap Stokes")
# ax3b.plot(subData25.ReP, subData25.scalDisConserv,
#           marker="o", ls='none', label="25 um Gap")
ax3b.set_xlabel('Reynolds Number')
ax3b.set_ylabel('Scalar Dissipation Rate (mol2/m3/s)')
ax3b.set_yscale('log')
ax3c.plot(subData100.ReP,
          subData100.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, color='k', label="Sim 100um")
ax3c.plot(subData25.ReP,
          subData25.dCdtMean/max(subData25.dCdtMean), ls='--',
          marker=None, color='k', label="Sim 25um")
ax3c.set_ylabel('Mean Reaction Rate Normalized to Max (.)')
ax3b.legend()

# Bonus figure: All trends superimposed
f10, ax10 = plt.subplots(ncols=1,nrows=1)
plt.subplots_adjust(right=0.75)
ax10.plot(subData100.ReP,
          subData100.dCdtMean/max(subData100.dCdtMean), ls='-',
          marker=None, label="Sim 100um")
ax10.plot(subData25.ReP,
          subData25.dCdtMean/max(subData25.dCdtMean), ls='-',
          marker=None, label="Sim 25um")
ax10.set_ylabel('Mean reaction rate normalized to max (.)')
ax10.set_xlabel('Reynolds Number')
ax10b = ax10.twinx()
ax10b.plot(subData100.ReP, subData100.MRT, marker="<",
           ls='none', label="100 um Gap")
ax10b.plot(subData25.ReP, subData25.MRT, marker="<",
           ls='none', label="25 um Gap")
ax10b.set_ylabel('Mean residence time (s)')
ax10b.set_yscale('log')
ax10c = ax10.twinx()
ax10c.spines.right.set_position(("axes", 1.2))

ax10c.plot(subData100.ReP, subData100.scalDisConserv,
          marker="o", ls='none', label="100 um Gap")
ax10c.plot(subData25.ReP, subData25.scalDisConserv,
          marker="o", ls='none', label="25 um Gap")
ax10c.set_yscale('log')
ax10c.set_ylabel('Scalar dissipation rate (mol2/m3/s)')

sns.despine(f2)
sns.despine(f3)

plt.show()