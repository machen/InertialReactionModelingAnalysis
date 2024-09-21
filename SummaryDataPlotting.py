import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#TODO: Update to pull different items for plotting

dataFile = "..\\SimulationSummaryData_v6.xlsx"
data = pd.read_excel(dataFile, sheet_name='Data',usecols="A:P,Z:AA,AK,AP,AR,BA,BE,BL,BM,BO,BS,BU",
                     names=["FileName", "r1", "r2", "d", "k",
                            "c", "Re", "ReP", "FlowCond", "uInlet", "EstRT",
                            "MRT", "MRTerr", "TimeRatio",
                            "PeNaive", "PePoreThroat",
                            "dCdtMean", "dCdtStd", "cLim","cTCPOLim", "cProdLim",
                            "DaAdvPoreThroat", "DaDiffPoreThroat",
                            "dCdtSum",
                            "dilConsLim","reacConsvLim",
                            "ReHd",
                            "ReDarcy"
                            ], skiprows=2)
data.sort_values(by=["FlowCond","ReP"], inplace=True)

plt.ion()
sns.set_context('poster', font_scale=1.25)
plt.rcParams['font.family'] = 'Cambria'
plt.rcParams['svg.fonttype'] = 'none'

# Color palettes. 1 is full color, 2 is muted
pal1 = sns.color_palette(as_cmap=True)
pal2 = sns.color_palette('muted', as_cmap=True)

subData100 = data.loc[(data.d == 100) & (data.k == 2000) &
                      (data.FlowCond == 'NS'), :]
subData25 = data.loc[(data.d == 25) & (data.k == 2000) &
                     (data.FlowCond == 'NS'), :]
subData50 = data.loc[(data.d == 50) & (data.k == 2000) &
                      (data.FlowCond == 'NS'), :]
# For Figure 2: Normalized Sum dCdt results, NS

rePlot = 'ReP'

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
ax1.plot(subData100[rePlot],
         subData100.dCdtSum/max(subData100.dCdtSum), ls='-',
         marker='o', label="Sim 100um NS", color=pal1[0])
ax1.plot(subData50[rePlot],
         subData50.dCdtSum/max(subData50.dCdtSum), ls='-',
         marker='o', label="Sim 50um NS", color=pal1[2])
ax1.plot(subData25[rePlot],
         subData25.dCdtSum/max(subData25.dCdtSum), ls='-',
         marker='o', label="Sim 25um NS", color=pal1[1])
# # Plot stokes data
# subData100Stokes = data.loc[(data.d == 100) & (data.k==2000) & (data.FlowCond=='Stokes'), :]
# ax1.plot(subData100Stokes[rePlot],
#          subData100Stokes.dCdtSum/max(subData100.dCdtSum), ls='--',
#          marker='o', label="Sim 100um Stokes")
# subData25Stokes = data.loc[(data.d == 25) & (data.k==2000) & (data.FlowCond=='Stokes'), :]
# ax1.plot(subData25Stokes[rePlot],
#          subData25Stokes.dCdtSum/max(subData25.dCdtSum), ls='--',
#          marker='o', label="Sim 25um Stokes")

ax1.legend()
ax1.set_xlabel('Reynolds Number')
ax1.set_ylabel('sum(dCdt) Normalized to Max (.)')

# ax1.set_xlim([-1, 105])
# ax1.set_ylim([0.3, 1.05])
# ax1.set_title('Check normalization values')
sns.despine(f1)

# For figure 3: Plots in mean tracer
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 20))
ax2.plot(subData100[rePlot], subData100.cLim, ls='-',
         marker='o', label="C_lim Tracer 100 um",
         color=pal1[0])
ax2.plot(subData100[rePlot], subData100.cTCPOLim, ls='--',
         marker='s', label="Mean TCPO 100 um",
         color=pal2[0])
ax2.plot(subData50[rePlot], subData50.cLim, ls='-',
         marker='o', label="C_lim Tracer 50 um",
         color=pal1[2])
ax2.plot(subData50[rePlot], subData50.cTCPOLim, ls='--',
         marker='s', label="Mean TCPO 50 um",
         color=pal2[2])
ax2.plot(subData25[rePlot], subData25.cLim, ls='-',
         marker='o', label="C_lim Tracer 25 um",
         color=pal1[1])
ax2.plot(subData25[rePlot], subData25.cTCPOLim, ls='--',
         marker='s', label="Mean TCPO 25 um",
         color=pal2[1])

ax2.legend()
ax2.set_xlabel('Reynolds Number')
ax2.set_ylabel('Concentration (mM)')
ax2.set_ylim(ymin=0.0, ymax=0.6)
sns.despine(f2)


# Figure 3: MRT
f3, ax3 = plt.subplots(ncols=1, nrows=1, figsize=(12, 10))
tReact = 1/3E-3/2000  # Half life assuming well mixed and equal reactants
# ax3c = ax3b.twinx()
ax3.plot(subData100[rePlot], subData100.MRT, marker="o",
         ls='none', label="100 um Gap", color=pal1[0])
ax3.plot(subData50[rePlot], subData50.MRT, marker="o",
         ls='none', label="50 um Gap", color=pal1[2])
ax3.plot(subData25[rePlot], subData25.MRT, marker="o",
         ls='none', label="25 um Gap", color=pal1[1])
ax3.set_xlabel('Reynolds Number')
ax3.set_ylabel('Mean residence time (s)')
ax3.set_yscale('log')
#ax3.set_ylim([0.0005, 500])
ax3.legend()
sns.despine(f3)

# Figure of Pe and relevant Da numbers
f4, ax4 = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

# Create data sets that give either the adv or diff Da depending on the Pe number
# TODO: This relies on Pe increasing monotonically from diffusive to advective and is really hacky. Should work aribitrarily. Maybe I should be actually just creating a "net Da" column
daAdv100 = subData100.loc[subData100.PePoreThroat > 1,'DaAdvPoreThroat']
daDiff100 = subData100.loc[subData100.PePoreThroat < 1, 'DaDiffPoreThroat']
da100 = pd.concat([daDiff100, daAdv100])
da100.sort_values(axis='index')

daAdv25 = subData25.loc[subData25.PePoreThroat > 1, 'DaAdvPoreThroat']
daDiff25 = subData25.loc[subData25.PePoreThroat < 1, 'DaDiffPoreThroat']
da25 = pd.concat([daDiff25, daAdv25])
da25.sort_values(axis='index')

daAdv50 = subData50.loc[subData50.PePoreThroat > 1, 'DaAdvPoreThroat']
daDiff50 = subData50.loc[subData50.PePoreThroat < 1, 'DaDiffPoreThroat']
da50 = pd.concat([daDiff50, daAdv50])
da50.sort_values(axis='index')

color1 = pal1[9]
color1desat = pal2[9]
color2 = pal1[3]
color2desat = pal2[3]
color3 = pal1[5]
color3desat = pal2[5]

ax4[0].plot(subData100[rePlot],
            subData100.PePoreThroat, ls='-',
            marker='o', label="Pe 100um",
            color=color1)
ax4[0].plot(subData50[rePlot],
            subData50.PePoreThroat, ls='-.',
            marker='^', label="Pe 50um",
            color=color3)
ax4[0].plot(subData25[rePlot],
            subData25.PePoreThroat, ls='--',
            marker='s', label="Pe 25um",
            color=color1desat)
ax4[0].plot([subData25[rePlot].min(), subData25[rePlot].max()],
            [1, 1], color='k', ls=':')
ax4[1].plot(subData100[rePlot],
            da100, ls = '-', color=color2,
            marker='o', label="Da 100 um")
ax4[1].plot(subData50[rePlot],
            da50, ls='-.', color=color3,
            marker='^', label="Da 50 um")
ax4[1].plot(subData25[rePlot],
            da25, ls='--', color=color2desat,
            marker='s', label="Da 25 um")
ax4[1].plot([subData25[rePlot].min(), subData25[rePlot].max()],
            [1, 1], color='k', ls=':')

# Also set axis properties
# ax4[0].spines['left'].set_color(color1)
# ax4[1].spines['left'].set_color(color1)
# # ax4[0].tick_params(axis='y', colors=color1, which='both')
# ax4[0].yaxis.label.set_color(color1)

# ax4[1].spines['right'].set_color(color2)
# # ax4[1].tick_params(axis='y',colors=color2, which='both')
# # ax4[1].yaxis.label.set_color(color2)

ax4[0].set_xlabel('Reynolds Number')
ax4[1].set_xlabel('Reynolds Number')
ax4[0].set_ylabel('Peclet Number (-)')
ax4[0].set_yscale('log')
ax4[1].set_ylabel('Dahmkohler Number (-)')
ax4[1].set_yscale('log')
ax4[0].legend()
ax4[1].legend()
ax4[0].set_title('Pore throat ND Numbers')
sns.despine(f4, right=False)


# Dilution Index
f5, ax5 = plt.subplots(1, 1,figsize=(12,10))
ax5.plot(subData100[rePlot], subData100.reacConsvLim, ls='-',
         marker='o', label="100 um",
         color=pal1[0])
ax5.plot(subData50[rePlot], subData50.reacConsvLim, ls='-',
         marker='o', label="50 um",
         color=pal1[2])
ax5.plot(subData25[rePlot], subData25.reacConsvLim, ls='-',
         marker='o', label="25 um",
         color=pal1[1])
ax5.set_xlabel('Reynolds Number')
ax5.set_ylabel('Reactor Ratio')
sns.despine(f5)

# Figure: Product in region where TCPO has crossed the boundary

f6, ax6 = plt.subplots(1, 1, figsize=(12,10))
ax6.plot(subData100[rePlot], subData100.cTCPOLim, ls='-', marker='o', label='100 um', color=pal1[0])
ax6.plot(subData50[rePlot], subData50.cTCPOLim, ls='-', marker='o', label='50 um', color=pal1[2])
ax6.plot(subData25[rePlot], subData25.cTCPOLim, ls='-',
         marker='o', label="25 um",
         color=pal1[1])
ax6.set_xlabel('Reynolds Number')
ax6.set_ylabel('Product (mM)')
ax6.set_title('Region where TCPO has crossed reactant boundary')
ax6.legend()

f7, ax7 = plt.subplots(1, 1, figsize=(12,10))
ax7.plot(subData100[rePlot], subData100.cProdLim/subData100.cLim,
         ls='-', marker='o', label='100 um', color=pal1[0])
# ax7.plot(subData50[rePlot], subData50.cProdLim/subData100.cLim,
#          ls='-', marker='o', label='50 um', color=pal1[2])
ax7.plot(subData25[rePlot], subData25.cProdLim/subData25.cLim,
         ls='-', marker='o', label="25 um",
         color=pal1[1])
ax7.set_xlabel('Reynolds Number')
ax7.set_ylabel('Availability Normalized Reaction Completion (-)')
ax7.set_title('Region where TCPO has crossed reactant boundary')
sns.despine(f7)
plt.show()