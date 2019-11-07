import pytools as pt
import numpy as np
from multiprocessing import Pool
import scipy.io

path_bulk = '/scratch/project_2000203/2D/BFH/reverted_ionosphere_field_boundary/'
path_save = '/users/dubartma/analysator/Data_Analysis/BFH/Data/'

bulkStart = 400
bulkEnd   = 1000
numproc   = 40

CellID = np.load(path_save+'CellID_SC.npy')

def get_field(step):
    
    bulkname = 'bulk.'+str(step).rjust(7,'0')+'.vlsv' 
    f        = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    xsize  = f.read_parameter('xcells_ini')
    xmax   = f.read_parameter("xmax")
    xmin   = f.read_parameter("xmin")
    dx     = (xmax-xmin)/xsize
    coords = f.get_cell_coordinates(CellID)
    coords = coords/dx

    B = f.read_fsgrid_variable('fg_b')
    B = B[int(coords[0]-0.5),0,int(coords[2]-0.5),:]

    return B    

pool         = Pool(numproc)
return_field = pool.map(get_field, range(bulkStart,bulkEnd+1))
return_field = np.array(return_field)

scipy.io.savemat(path_save+'BFH.mat',dict(B=return_field))
print(path_save+'BFH.mat')
