# Determine angle between z-axis and B-field for FFT

import pytools as pt
import numpy as np
from multiprocessing import Pool

def get_angle(step)

    f = pt.vlsvfile.VlsvReader(

path_bulk = '/scratch/project_2000203/2D/BFH/reverted_ionosphere_field_boundary/'
path_save = '/users/dubartma/BFH/Data/'

bulkStart = 200
bulkEnd   = 400
bulkTot   = bulkEnd - bulkStart +1 
RE        = 6371e3 # m
numproc   = 80

bulkname  = 'bulk.'+str(bulkStart).rjust(7,'0')+'.vlsv'
print(bulkname)

f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    # Get size of simulation (in cells)
xsize = f.read_parameter('xcells_ini')
xmax = f.read_parameter("xmax")
xmin = f.read_parameter("xmin")
dx   = (xmax-xmin)/xsize

zsize = f.read_parameter('zcells_ini')
if zsize == 1:
    plane = 'equatorial'
    ysize = f.read_parameter('ycells_ini')
    ymax = f.read_parameter("ymax")
    ymin = f.read_parameter("ymin")
    dy   = (ymax-ymin)/ysize
else:
    plane = 'polar'
    ysize = zsize
    ymax = f.read_parameter("zmax")
    ymin = f.read_parameter("zmin")
    dy   = (ymax-ymin)/ysize
print('-----------------------------')
print('Run: '+run+', plane: '+plane)
print('-----------------------------')

    # Box coordinates # m
x1 = int((3.0 *RE - xmin) / dx)
x2 = int((6.0 *RE - xmin) / dx)
y1 = int((15.0 *RE - ymin) / dy)
y2 = int((18.0 *RE - ymin) / dy)

# One point in the middle
xmid = int((x2 + x1) / 2)
ymid = int((y2 + y1) / 2)

# Initialise field
Bmid = np.zeros([3,bulkTot])

print('Start processing files')
# Computes angle between magnetic field and z-axis (or y-axis for equatorial plane)
# ASSUME THAT THE ANGLE DOESNT VARY TOO MUCH DURING THE ANALYSIS
if __name__ == 'main':
    pool         = Pool(numproc)
    return_angle = pool.map(get_angle, range(bulkStart,bulkStop+1))
else:
    print("Problem with angle parallelization")
print('All files processed')


# Time avg of B per component & norm
B_avg  = np.mean(Bmid,axis=1) 
B_norm = np.linalg.norm(B_avg)

# Angle
ez = [0,0,1]
angle = abs(np.arccos(np.dot(B_avg,ez)/B_norm) * 180.0 / np.pi -180.0) # Degrees

np.save(path_save+'angle.npy',angle)
print('Angle = '+str(angle)+', saved.')
