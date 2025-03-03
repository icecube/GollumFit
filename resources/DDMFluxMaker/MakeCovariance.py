import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.interpolate import fitpack2
sys.modules['scipy.interpolate._fitpack2'] = sys.modules['scipy.interpolate.fitpack2']

def _remap_covariance(original_order, new_order, cov):
    cov_new = np.zeros((len(new_order), len(new_order)))
    remap = original_order
    for i in range(cov_new.shape[0]):
        for j in range(cov_new.shape[1]):
            cov_new[i, j] = cov[remap[new_order[i]], remap[new_order[j]]]
    return cov_new


syst_dict = {
'K+_158G'  : 'he_K+',
'K-_158G'  : 'he_K-',
'pi+_20T'  : 'vhe1_pi+',
'pi-_20T'  : 'vhe1_pi-',
'K+_2P'    : 'vhe3_K+',
'K-_2P'    : 'vhe3_K-',
'pi+_2P'   : 'vhe3_pi+',
'pi-_2P'   : 'vhe3_pi-',
'p_2P'     : 'vhe3_p',
'n_2P'     : 'vhe3_n',
'GSF_1'    : 'GSF_1',
'GSF_2'    : 'GSF_2',
'GSF_3'    : 'GSF_3',
'GSF_4'    : 'GSF_4',
'GSF_5'    : 'GSF_5',
'GSF_6'    : 'GSF_6'}


flux = pickle.load(open('./daemonflux/daemonsplines_IceCube_202303_1.pkl','rb'),encoding="latin1")
cal = pickle.load(open('./daemonflux/daemonsplines_calibration_202303_1.pkl','rb'),encoding="latin1")

known_params = flux[0]

print(known_params)

# exclude_list = ["pbar", "pi+", "pi-", "K+", "K-", "N"]
# if "he_p" in known_params: exclude_list += ["p", "n"]
# for p in known_params:
#     if "pbar" in p: exclude_list += [p]
#     if "le_" in p: exclude_list += [p]
# exclude_list += ['he_pi+','he_pi-','he_p','he_n','vhe3_p','vhe3_n']
# known_params = [k for k in known_params if k not in exclude_list]

known_params = list(syst_dict.keys())


original_param_order = {}
for ip, p in enumerate(cal["cov_params"]) :
    if p in known_params :  
        original_param_order[p] = ip

print('original_param_order',original_param_order)
print('known_params',known_params)

n_physics_params = max(original_param_order.values()) + 1

# Create a new covariance with the correct order of parameters
cal_corr = _remap_covariance(original_param_order, known_params, cal["corr_matrix"][:n_physics_params, :n_physics_params])

np.save("./flux_correlations",cal_corr)

fig = plt.figure(figsize=(12,8.15))
ax1 = plt.subplot()
img = plt.imshow(cal_corr, interpolation='nearest', cmap='RdBu')
plt.clim(-1,1)
cbar = fig.colorbar(img, ax=ax1)
cbar.ax.set_ylabel('Corr')
plt.xticks(range(len(syst_dict)),syst_dict.values(),rotation=90)
plt.yticks(range(len(syst_dict)),syst_dict.values())
plt.tight_layout()
plt.savefig("./ddm_corr.pdf")
plt.close()
