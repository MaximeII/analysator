# Save arrays to text file

import numpy as np
import scipy.io

path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/V2/BCE/Data/'

B   = np.load(path_save+'B_SC.npy')
E   = np.load(path_save+'E_SC.npy')
v   = np.load(path_save+"v_SC.npy")
rho = np.load(path_save+"rho_SC.npy")
#T_SC        = np.load(path_save+"T_SC.npy")
#Taniso_SC   = np.load(path_save+'Taniso_SC.npy')
#Tpara_SC    = np.load(path_save+"Tpara_SC.npy")
#Tperp_SC    = np.load(path_save+"Tperp_SC.npy")
#betaPara_SC = np.load(path_save+"betaPara_SC.npy")
#betaPerp_SC = np.load(path_save+"betaPerp_SC.npy")
#J           = np.load(path_save+'J.npy')
#Jperp       = np.load(path_save+'Jperp.npy')
#vperp       = np.load(path_save+'vperp.npy')
#vpara       = np.load(path_save+'vpara.npy')

scipy.io.savemat(path_save+'BCE.mat',dict(B=B,E=E,v=v,rho=rho))

#scipy.io.savemat(path_save+'Mojtaba.mat',dict(B=B,E=E,v=v_SC,n=rho_SC,T=T_SC,Taniso=Taniso_SC,Tpara=Tpara_SC,Tperp=Tperp_SC,betaPara=betaPara_SC,betaPerp=betaPerp_SC,J=J,Jperp=Jperp,vperp=vperp,vpara=vpara))
