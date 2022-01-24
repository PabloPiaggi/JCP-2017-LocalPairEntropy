import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm2
from mpl_toolkits.basemap import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate

###################################################################
# Plot options
###################################################################
font = {'family' : 'serif',
        'serif'   : 'palatino',
        # 'sans-serif'    : 'Computer Modern Sans serif', 
        'style'   : 'normal',
        'variant'   : 'normal',
        'stretch'   : 'normal',
        'weight'   : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = (8.3, 5.8)

###################################################################
# Set constants
###################################################################

###################################################################
# Colors
###################################################################
# Palette 1
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
MarkerSize=20
palette=tableau20

###################################################################
# Plot data
###################################################################

#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()

#freeEnergy=np.genfromtxt("../../results2.txt")
#freeEnergy[:,2] += 0.032
#committor=np.genfromtxt("../results.txt")
#coordinate = -np.linspace(0,freeEnergy.shape[0],freeEnergy.shape[0])

#ax2.plot([-12.57,-12.57],[-0.5,1.5],'-',color="grey",linewidth=2.0,alpha=0.8)
#ax1.plot(coordinate,freeEnergy[:,2],'-',color=palette[0],linewidth=4.0,alpha=0.8)
#ax1.fill_between(coordinate,0,freeEnergy[:,2],color=palette[0],alpha=0.1)
#ax2.plot(coordinate,committor[:,2],'--',color=palette[2],linewidth=4.0,alpha=0.8)

liquid1=np.genfromtxt('liquid/histoS')
fcc1=np.genfromtxt('fcc/histoS')
liquid2=np.genfromtxt('liquid/histoSa')
fcc2=np.genfromtxt('fcc/histoSa')

plt.plot(liquid1[:,0],liquid1[:,1],'--',color=palette[2],linewidth=2.0,alpha=0.6)
plt.fill_between(liquid1[:,0],0,liquid1[:,1],color=palette[2],alpha=0.1)
plt.plot(fcc1[:,0],fcc1[:,1],'--',color=palette[0],linewidth=2.0,alpha=0.6)
plt.fill_between(fcc1[:,0],0,fcc1[:,1],color=palette[0],alpha=0.1)
plt.plot(liquid2[:,0],liquid2[:,1],color=palette[2],linewidth=2.0,alpha=0.8,label="liquid")
plt.fill_between(liquid2[:,0],0,liquid2[:,1],color=palette[2],alpha=0.2)
plt.plot(fcc2[:,0],fcc2[:,1],color=palette[0],linewidth=2.0,alpha=0.8,label="fcc")
plt.fill_between(fcc2[:,0],0,fcc2[:,1],color=palette[0],alpha=0.2)

# Now add the legend with some customizations.
legend = plt.legend(loc='upper left')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(18)

###################################################################
# Plot options
###################################################################

plt.xlim([-5,-1.])
#plt.ylim([0,3.5])
#ax1.set_ylim([-0.2,22])
#ax2.set_ylim([-0.01,1.01])
plt.xlabel('Pair entropy (k$_B$)')
plt.ylabel('Probability density')
#ax1.set_xlabel('Minimum free energy path')
#ax2.set_xlabel('Minimum free energy path')
#ax1.set_ylabel('Free energy (kJ/mol)')
#ax2.set_ylabel('Committor $p_s$')
plt.tick_params(axis='x', pad=10)
plt.tick_params(axis='y', pad=10)
#plt.xticks([-30,0],["Liquid","Solid"]) 
#plt.xticks(np.linspace(-7,-1,4),np.linspace(-7,-1,4).astype(int))
#plt.yticks(np.linspace(0,4,5),np.linspace(0,4,5).astype(int))
#plt.show()
plt.savefig('histograms2.pdf', bbox_inches='tight')
plt.savefig('histograms2.png', bbox_inches='tight', dpi=300)

