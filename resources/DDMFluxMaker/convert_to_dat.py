import numpy as np
import matplotlib.pyplot as plt
from daemonflux import Flux
import pickle

import sys
from scipy import interpolate
from scipy.interpolate import fitpack2
sys.modules['scipy.interpolate._fitpack2'] = sys.modules['scipy.interpolate.fitpack2']



# Making a 2D flux map (E vs zenith) for a particle flux
def make2Dmap(fname, particle, vegrid, vparams={}):
    
    # Point to the files
    quickflux = Flux(spl_file='./daemonflux/daemonsplines_IceCube_202303_1.pkl', cal_file='./daemonflux/daemonsplines_calibration_202303_1.pkl', use_calibration=True)

    icangles = list(quickflux.zenith_angles)
    
    # Sorting the angles in the dictionary
    icangles_array = np.array(icangles, dtype=float)
    mysort = icangles_array.argsort()
    icangles = np.array(icangles)[mysort][::-1]

    flux_ref = np.zeros([len(vegrid), len(icangles)])

    costheta_angles = np.zeros(len(icangles))

    for index in range(len(icangles)):
        costheta = np.cos(np.deg2rad(np.float(icangles[index])))
        costheta_angles[index] = costheta
        if   fname=='conv' : flux_ref[:,index] = quickflux.flux(vegrid, icangles[index], particle, params=vparams)
        elif fname=='pr'   : flux_ref[:,index] = quickflux.flux(vegrid, icangles[index], 'total_'+particle, params=vparams) - quickflux.flux(vegrid, icangles[index], particle, params=vparams)

    spl = interpolate.RectBivariateSpline(np.log10(vegrid), costheta_angles, flux_ref)  

    return spl


def makefile(fname,sname,sout) :

    params_dict = {}
    if sname!='bestfit' : params_dict[sname] = 1

    print(params_dict)

    ne = 350
    nc = 100

    energies = np.logspace(2,6,ne)

    spl_flux = {}
    for pname in ['numu','antinumu','nue','antinue'] : spl_flux[pname] = make2Dmap(fname=fname, particle=pname, vegrid=energies, vparams=params_dict)

    cos_zeniths = np.linspace(-1, .2, nc)

    out_flux = np.empty([ne*nc,8])
                
    e_count = 0
    c_count = 0
    for i in range(ne*nc):
        energy = energies[e_count]
        cosz   = cos_zeniths[c_count]

        conv_numu    = spl_flux['numu'](np.log10(energy), cosz) / np.power(energy,3)
        conv_anumu   = spl_flux['antinumu'](np.log10(energy), cosz) / np.power(energy,3)
        conv_nue     = spl_flux['nue'](np.log10(energy), cosz) / np.power(energy,3)
        conv_anue    = spl_flux['antinue'](np.log10(energy), cosz) / np.power(energy,3)
        conv_nutau   = 0.
        conv_anutau  = 0.
        
        out_flux[i, 0] = cosz
        out_flux[i, 1] = energy
        out_flux[i, 2] = conv_nue
        out_flux[i, 3] = conv_anue
        out_flux[i, 4] = conv_numu
        out_flux[i, 5] = conv_anumu
        out_flux[i, 6] = conv_nutau 
        out_flux[i, 7] = conv_anutau

        e_count += 1

        if e_count % 350 == 0:
            c_count += 1
            e_count = 0

    if sname=='bestfit' : np.savetxt('./dat_files/v0.6.0/ddm_'+fname+'_'+sout+'.dat', out_flux)
    else                : np.savetxt('./dat_files/v0.6.0/ddm_'+sout+'.dat', out_flux)


syst_dict = {
#le
'K+_31G'   : 'le_K+' ,
'K-_31G'   : 'le_K-' ,
'pi+_31G'  : 'le_pi+' ,
'pi-_31G'  : 'le_pi-' ,
#he
'K+_158G'  : 'he_K+',
'K-_158G'  : 'he_K-',
'pi+_158G' : 'he_pi+',
'pi-_158G' : 'he_pi-',
'n_158G'   : 'he_n',
'p_158G'   : 'he_p',
#vhe1
'pi+_20T'  : 'vhe1_pi+',
'pi-_20T'  : 'vhe1_pi-',
#vhe3
'K+_2P'    : 'vhe3_K+',
'K-_2P'    : 'vhe3_K-',
'n_2P'     : 'vhe3_n',
'p_2P'     : 'vhe3_p',
'pi+_2P'   : 'vhe3_pi+',
'pi-_2P'   : 'vhe3_pi-',
#gsf
'GSF_1' : 'GSF_1',
'GSF_2' : 'GSF_2',
'GSF_3' : 'GSF_3',
'GSF_4' : 'GSF_4',
'GSF_5' : 'GSF_5',
'GSF_6' : 'GSF_6'}


makefile('conv','bestfit','bestfit')
makefile('pr','bestfit','bestfit')

for sname,sout in syst_dict.items() : makefile('conv',sname,sout)
