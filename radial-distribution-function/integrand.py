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

def calcIntegrand(fcc):
	integrand=np.zeros(fcc[:,0].shape[0])
	for i in range(0,fcc[:,0].shape[0],1):
		if (fcc[i,1]<1.e-10):
			integrand[i]=fcc[i,0]*fcc[i,0]
		else:
			integrand[i]=(fcc[i,1]*np.log(fcc[i,1])-fcc[i,1] +1)*fcc[i,0]*fcc[i,0]
	return integrand

###################################################################
# Plot data
###################################################################

f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)

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

fcc=np.genfromtxt('Fcc/Gofr/gofr.dat')
hcp=np.genfromtxt('Hcp/Gofr/gofr.dat')
bcc=np.genfromtxt('Bcc/Gofr/gofr.dat')
liquid=np.genfromtxt('Liquid/Gofr/gofr.dat')

ax1.plot(fcc[:,0]/10.,calcIntegrand(fcc),color=palette[4],linewidth=2.0,alpha=0.8,label="fcc")
ax1.fill_between(fcc[:,0]/10.,0.,calcIntegrand(fcc),color=palette[4],linewidth=2.0,alpha=0.15)
ax1.text(1.55,10.,"fcc", rotation='vertical')
ax2.plot(hcp[:,0]/10.,calcIntegrand(hcp),color=palette[6],linewidth=2.0,alpha=0.8,label="hcp")
ax2.fill_between(hcp[:,0]/10.,0.,calcIntegrand(hcp),color=palette[6],linewidth=2.0,alpha=0.15)
ax2.text(1.55,10.,"hcp", rotation='vertical')
ax3.plot(bcc[:,0]/10.,calcIntegrand(bcc),color=palette[0],linewidth=2.0,alpha=0.8,label="bcc")
ax3.fill_between(bcc[:,0]/10.,0.,calcIntegrand(bcc),color=palette[0],linewidth=2.0,alpha=0.15)
ax3.text(1.55,10.,"bcc", rotation='vertical')
ax4.plot(liquid[:,0]/10.,calcIntegrand(liquid),color=palette[2],linewidth=2.0,alpha=0.8,label="liquid")
ax4.fill_between(liquid[:,0]/10.,0.,calcIntegrand(liquid),color=palette[2],linewidth=2.0,alpha=0.15)
ax4.text(1.55,13.,"liquid", rotation='vertical')

#plt.plot(bcc2[:,0],bcc2[:,1],color=palette[0],linewidth=2.0,alpha=0.8,label="bcc")
#plt.fill_between(bcc2[:,0],0,bcc2[:,1],color=palette[0],alpha=0.2)

# Now add the legend with some customizations.
#legend = plt.legend(loc='upper left')

# Set the fontsize
#for label in legend.get_texts():
#    label.set_fontsize(18)

###################################################################
# Plot options
###################################################################

plt.xlim([0,1.5])
#plt.ylim([0,3.5])
#ax1.set_ylim([-0.2,22])
#ax2.set_ylim([-0.01,1.01])
plt.xlabel('Distance (nm)')
#plt.ylabel('Probability density')
#ax2.set_ylabel('$\{ g(r)\log[g(r)]-g(r)+1 \}r^2$')
ax3.text(-0.15, 22., '$\{ g(r)\log[g(r)]-g(r)+1 \}r^2$', ha='center', va='center', rotation='vertical')
#ax2.set_xlabel('Minimum free energy path')
#ax1.set_ylabel('Free energy (kJ/mol)')
#ax2.set_ylabel('Committor $p_s$')
plt.tick_params(axis='x', pad=10)
plt.tick_params(axis='y', pad=10)
ax1.set_yticks([0,10,20]) 
ax1.set_yticklabels(np.array([0,10,20]).astype(int)) 
ax2.set_yticks([0,10,20]) 
ax2.set_yticklabels(np.array([0,10,20]).astype(int)) 
ax3.set_yticks([0,10,20]) 
ax3.set_yticklabels(np.array([0,10,20]).astype(int)) 
ax4.set_yticks([0,10,20]) 
ax4.set_yticklabels(np.array([0,10,20]).astype(int)) 
#plt.xticks(np.linspace(-7,-1,4),np.linspace(-7,-1,4).astype(int))
#plt.yticks(np.linspace(0,4,5),np.linspace(0,4,5).astype(int))
#plt.show()
plt.savefig('integrand.pdf', bbox_inches='tight')
plt.savefig('integrand.png', bbox_inches='tight', dpi=300)

