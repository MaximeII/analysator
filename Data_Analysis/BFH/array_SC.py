# Get data from 3 spacecraft

import pytools as pt
import numpy as np
import sys

path_bulk = '/proj/vlasov/2D/BCE/bulk/'
path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/V2/BCE/Data/'

# Load data
CellID = int(np.load(path_save+"CellID_SC.npy"))
print "Data loaded"

if len(sys.argv)==3: # Starting and end bulk given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting bulk given, generate one bulk
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    #Open file
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    B_SC        = f.read_variable('B',CellID)
    E_SC        = f.read_variable('E',CellID)
    rho_SC      = f.read_variable('rho',CellID)
    T_SC        = f.read_variable('Temperature',CellID)
    Tpara_SC    = f.read_variable('TParallel',CellID)
    Tperp_SC    = f.read_variable('TPerpendicular',CellID)
    Taniso_SC   = f.read_variable('TPerpOverPar',CellID)
    betaPara_SC = f.read_variable('betaParallel',CellID)
    betaPerp_SC = f.read_variable('betaPerpendicular',CellID) 
    v_SC        = np.array(f.read_variable('v',CellID))    

    # Save
    np.save(path_save+"B_SC_"+str(j)+".npy",B_SC)
    np.save(path_save+"E_SC_"+str(j)+".npy",E_SC)
    np.save(path_save+"rho_SC_"+str(j)+".npy",rho_SC)
    np.save(path_save+"T_SC_"+str(j)+".npy",T_SC)
    np.save(path_save+"Tpara_SC_"+str(j)+".npy",Tpara_SC)
    np.save(path_save+"Tperp_SC_"+str(j)+".npy",Tperp_SC)
    np.save(path_save+"Taniso_SC_"+str(j)+".npy",Taniso_SC)
    np.save(path_save+"betaPara_SC_"+str(j)+".npy",betaPara_SC)
    np.save(path_save+"betaPerp_SC_"+str(j)+".npy",betaPerp_SC)
    np.save(path_save+"v_SC_"+str(j)+".npy",v_SC)
