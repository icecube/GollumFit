import numpy as np
import matplotlib.pyplot as plt
import sys

def print_matrix(mname,m) :
    print(mname+' :')
    for i in range(len(m)) : 
        txt = ''
        for j in range(len(m)) : txt += "%+0.5f " % m[i,j]
        print(txt)

def correlation(cov) :
    corr = np.zeros_like(cov)
    for i in range(len(cov)) : 
        for j in range(len(cov)) : 
            corr[i,j] = cov[i,j]/np.sqrt(cov[i,i])/np.sqrt(cov[j,j])
    return corr


inwidths = sys.argv[1]
outdir   = sys.argv[2]


Modes=5

ampamp = np.loadtxt(inwidths+"Widths/AmpAmpWidths.txt")
phsphs = np.loadtxt(inwidths+"Widths/PhsPhsWidths.txt")
ampphs = np.loadtxt(inwidths+"Widths/AmpPhsWidths.txt")

Matrix=np.ndarray(shape=(Modes+Modes-1,Modes+Modes-1))
Matrix=np.zeros_like(Matrix)
for i in range(0, Modes):
    for j in range(0, i):
        mask1=(ampamp[:,0]==i)*(ampamp[:,1]==j)
        Matrix[i,j]=(ampamp[mask1][0][4]-ampamp[mask1][0][2])

for i in range(1, Modes):
    for j in range(1, i):
        mask1=(phsphs[:,0]==i)*(phsphs[:,1]==j)
        Matrix[i+Modes-1,j+Modes-1]=phsphs[mask1][0][4]-phsphs[mask1][0][2]


for i in range(1, Modes):
    for j in range(0, Modes):
        mask1=(ampphs[:,0]==j)*(ampphs[:,1]==i)
        Matrix[i+Modes-1,j]=ampphs[mask1][0][4]-ampphs[mask1][0][2]
        
for i in range(0,len(Matrix)):
    for j in range(0, i):
        Matrix[j,i]=Matrix[i,j]
    Matrix[i,i]=1


NuisHessian=np.ones_like(Matrix)
for i in range(0,2*Modes-1):
    for j in range(0,2*Modes-1):
        if (i != j ):
            NuisHessian[i,j]= Matrix[i,j]**-2-1.

CovNuis=np.linalg.inv(NuisHessian)


Amp1sig=np.loadtxt(inwidths+"Widths/AmpWidths.txt")
Phs1sig=np.loadtxt(inwidths+"Widths/PhsWidths.txt")
AmpWidths=Amp1sig[0:Modes,3]-Amp1sig[0:Modes,1]
PhsWidths=Phs1sig[0:Modes-1,3]-Phs1sig[0:Modes-1,1]
Widths = np.concatenate((AmpWidths,PhsWidths))

print_matrix('Covariance',CovNuis)

Realization = np.random.multivariate_normal(np.zeros(9),CovNuis,check_valid='ignore',size=100000000)
Realization = np.transpose(Realization)

RealCov = np.cov(Realization)
RealCorr = correlation(RealCov)

print_matrix('Covariance',RealCov)
print_matrix('Correlations',RealCorr)
print('det:',np.linalg.det(RealCorr))
print('eig:',np.linalg.eigvalsh(RealCorr))

# np.save(outdir+"/icegrad_correlations",RealCorr)

XTickNames=[]
for i in range(0,Modes) : XTickNames.append("Amp "+str(i))
for i in range(1,Modes) : XTickNames.append("Phs "+str(i))

fig = plt.figure(figsize=(12,8.15))
ax1 = plt.subplot()
img = plt.imshow(RealCorr, interpolation='nearest', cmap='RdBu')
plt.clim(-1,1)
plt.plot([-0.5,2.*Modes-1.5],[Modes-0.5,Modes-0.5],'--',color='black')
plt.plot([Modes-0.5,Modes-0.5],[-0.5,2.*Modes-1.5],'--',color='black')
XTickNames=[]
for i in range(0,Modes) : XTickNames.append("Amp "+str(i))
for i in range(1,Modes) : XTickNames.append("Phs "+str(i))
plt.xticks(range(2*Modes-1),XTickNames,rotation=90)
plt.yticks(range(2*Modes-1),XTickNames)
plt.xlim(-0.5,2*Modes-1.5)
plt.ylim(2*Modes-1.5,-0.5)
cbar = fig.colorbar(img, ax=ax1)
cbar.ax.set_ylabel('Ratio')
plt.tight_layout()
plt.savefig("/n/home05/agarciasoto/plots_golemfit/ice_corr.pdf")
plt.close()




