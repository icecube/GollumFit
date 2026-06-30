import h5py
import numpy as np
import sys
import os

sys.path.append(os.getenv('GOLEMSPACE')+"/sources/nuSQuIDS/lib")
import nuSQuIDS as nsq


def get_weight(file,w) :
    weight = np.zeros(len(file["NuEnergy"]['value']))
    nuenergy = np.array(file["NuEnergy"]['value'])*1e9 #in eV
    nucth    = np.cos(np.array(file["NuZenith"]['value']))
    nupdg    = np.array(file["PrimaryType"]['value'])
    onewgt   = np.array(file["oneweight"]['value'])
    flavor   = (np.abs(nupdg)-12)/2
    neutype  = nupdg<0
    for i in range(len(nuenergy)) : 
        weight[i] += w.EvalFlavor(int(flavor[i]),float(nucth[i]),float(nuenergy[i]),int(neutype[i]),False)*onewgt[i]
    return weight



topology  = sys.argv[1] # Start_X, Through_X, All
varname   = sys.argv[2] # MultisimAmplitudes, MultisimPhases
splitmode = int(sys.argv[3]) # 0-12
influxdir = sys.argv[4]
outdir    = sys.argv[5]

# READ FILES
allfilenames = {}
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/LE/Ares/Merged/Ares_IC86.AVG_1.27_BDT_6867.h5"] = 6867
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/19.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/18.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/17.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/16.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/15.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/14.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/13.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/12.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/11.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/10.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/9.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/8.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/7.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/6.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/5.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/4.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/3.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/2.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/1.h5"] = 19800
allfilenames["/n/holylfs05/LABS/arguelles_delgado_lab/Lab/meows_dnn/mc/Systematics/Multisim/HE/Ares/Merged/Ares_IC86.AVG_1.27_BDT/0.h5"] = 19800
nebins = 24
emin   = 300
emax   = 1e5


nusq_conv   = nsq.nuSQUIDSAtm(influxdir+'/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5')
nusq_prompt = nsq.nuSQUIDSAtm(influxdir+'/prompt_atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5')
nusq_astro  = nsq.nuSQUIDSAtm(influxdir+'/astro_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5')


ene = np.array([], dtype=float)
cth = np.array([], dtype=float)
wc  = np.array([], dtype=float)
wp  = np.array([], dtype=float)
wa  = np.array([], dtype=float)
cut = np.array([], dtype=bool)
var = np.array([], dtype=float)
for filename, nfiles in allfilenames.items() :
    file = h5py.File(filename,'r')

    print(filename)
    
    var = np.append(var,file[varname][:,splitmode])
    ene = np.append(ene,file['DnnEnergy']['value'])
    cth = np.append(cth,np.cos(file['MuExZenith']['value']))    
    wc  = np.append(wc,get_weight(file,nusq_conv)/nfiles)
    wp  = np.append(wp,get_weight(file,nusq_prompt)/nfiles)
    wa  = np.append(wa,get_weight(file,nusq_astro)*4.72/6./nfiles)
    c1  = np.cos(np.array(file['MuExZenith']['value']))<0.
    c2  = np.array(file['DnnEnergy']['value'])>emin
    c3  = np.array(file['DnnEnergy']['value'])<emax
    if   'Start' in topology   : c4 = np.array(file['DeepStart']['value'])>float(topology.split("_")[-1])
    elif 'Through' in topology : c4 = np.array(file['DeepStart']['value'])<float(topology.split("_")[-1])
    elif 'All' in topology     : c4 = np.cos(np.array(file['MuExZenith']['value']))<0.
    cut = np.append(cut,c1 & c2 & c3 & c4)



var = var[cut]
ene = ene[cut]
cth = cth[cut]
wc  = wc[cut]
wp  = wp[cut]
wa  = wa[cut]

print(np.sum(wc)*1e6,np.sum(wp)*1e6,np.sum(wa)*1e6)

wgt = wc + wp + wa

print("Splitting Mode:",str(splitmode))
       
cut_pos = var>0
cut_neg = var<0       

mc_ene_pos = ene[ cut_pos ]
mc_ene_neg = ene[ cut_neg ]
mc_cth_pos = cth[ cut_pos ]
mc_cth_neg = cth[ cut_neg ]
mc_wgt_pos = wgt[ cut_pos ]
mc_wgt_neg = wgt[ cut_neg ]

print(np.sum(wgt)*1e6,np.sum(mc_wgt_pos)*1e6,np.sum(mc_wgt_neg)*1e6)

# SET BINNING 
ene_bins = np.logspace(np.log10(emin),np.log10(emax),nebins+1)
cth_bins = np.linspace(-1,0.,21)


######## CALCULATE WEIGHTED UNCERTAINTIES

errors_pos, x_edges, y_edges = np.histogram2d( mc_cth_pos, mc_ene_pos, bins=[cth_bins,ene_bins], weights=np.power(mc_wgt_pos,2))
errors_neg, x_edges, y_edges = np.histogram2d( mc_cth_neg, mc_ene_neg, bins=[cth_bins,ene_bins], weights=np.power(mc_wgt_neg,2))

errors_pos  = np.sqrt(errors_pos)
errors_neg  = np.sqrt(errors_neg)


########### GET MC DISTRIBUTIONS ##########

dist_pos, x_edges, y_edges = np.histogram2d( mc_cth_pos, mc_ene_pos, bins=[cth_bins,ene_bins], weights=mc_wgt_pos)
dist_neg, x_edges, y_edges = np.histogram2d( mc_cth_neg, mc_ene_neg, bins=[cth_bins,ene_bins], weights=mc_wgt_neg)       

ene_bin_lowedge = ene_bins[:-1]
ene_bin_upedge  = ene_bins[1:]
cth_bin_lowedge = cth_bins[:-1]
cth_bin_upedge  = cth_bins[1:]

file = open(outdir+"/"+topology+"_SplitCounts_"+varname+"_"+str(splitmode)+".csv", "w")
for i in range(len(cth_bin_lowedge)) : 
    for j in range(len(ene_bin_lowedge)) : 
        file.write(str(cth_bin_lowedge[i])+" "+str(cth_bin_upedge[i])+" "+str(ene_bin_lowedge[j])+" "+str(ene_bin_upedge[j])+" "+str(dist_pos[i][j])+" "+str(errors_pos[i][j])+" "+str(dist_neg[i][j])+" "+str(errors_neg[i][j])+'\n')
file.close()

print("Done! :-) \n")






