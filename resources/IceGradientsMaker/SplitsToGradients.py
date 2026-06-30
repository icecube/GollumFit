import numpy as np
import sys

topology  = sys.argv[1] # Start, Through
insplits  = sys.argv[2]
infrac    = sys.argv[3]
inwidths  = sys.argv[4]
outgrad   = sys.argv[5]

print(topology)
print(insplits)
print(infrac)
print(inwidths)
print(outgrad)


def SortIt(a):
    a=np.transpose(a)
    return np.transpose(a[::, a[0,].argsort()[::]])

GenScaleAmp=np.loadtxt(infrac+"/GenerationScaleAmp.dat")
GenScalePhs=np.loadtxt(infrac+"/GenerationScalePhs.dat")

OneSigmasAmpUnsr=np.loadtxt(inwidths+"/AmpWidths.txt")
OneSigmasPhsUnsr=np.loadtxt(inwidths+"/PhsWidths.txt")
OneSigmasAmp=SortIt(OneSigmasAmpUnsr)
OneSigmasPhs=SortIt(OneSigmasPhsUnsr)

CorrFactorAmp=np.zeros(10)
CorrFactorPhs=np.zeros(10)
for i in range(0,10): CorrFactorAmp[i] = ((OneSigmasAmp[i,3]-OneSigmasAmp[i,2])/2.)/GenScaleAmp[1,i]*np.sqrt(3.142/2.)
for i in range(1,10): CorrFactorPhs[i] = ((OneSigmasPhs[i-1,3]-OneSigmasPhs[i-1,2])/2.)/GenScalePhs[1,i]*np.sqrt(3.142/2.)
   
print("CorrFactorAmp",CorrFactorAmp)
print("CorrFactorPhs",CorrFactorPhs)

PhsGrads=[]
for grad in range(1,10):
    Phs=np.loadtxt(insplits+"/"+topology+"_SplitCounts_MultisimPhases_"+str(grad)+".csv")
    for line in Phs:
        if((line[4]+line[2])>0):
            PhsGrads.append([grad,line[0],line[1],line[2],line[3],(line[4]-line[6])/(line[6]+line[4])*CorrFactorPhs[grad],(line[7]**2+line[5]**2)**0.5/(line[6]+line[4])*CorrFactorPhs[grad]])


AmpGrads=[]
for grad in range(0,10):
    Amp=np.loadtxt(insplits+"/"+topology+"_SplitCounts_MultisimAmplitudes_"+str(grad)+".csv")
    for line in Amp:
        if((line[4]+line[2])>0):
            AmpGrads.append([grad,line[0],line[1],line[2],line[3],(line[4]-line[6])/(line[6]+line[4])*CorrFactorAmp[grad],(line[7]**2+line[5]**2)**0.5/(line[6]+line[4])*CorrFactorAmp[grad]])

np.savetxt(outgrad+"/"+topology+"_PhsGrads.csv",PhsGrads)
np.savetxt(outgrad+"/"+topology+"_AmpGrads.csv",AmpGrads)


