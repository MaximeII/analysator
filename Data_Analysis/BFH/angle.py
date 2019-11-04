# Determine angle between z-axis and B-field for FFT

import pytools as pt
import numpy as np

path_bulk = '/scratch/project_2000203/2D/BFH/reverted_ionosphere_field_boundary/'
path_save = '/users/dubartma/BFH/Data/'

bulkStart = 200
bulkEnd   = 400
bulkTot   = bulkEnd - bulkStart +1 

# Box coordinates # m 
bulkname  = 'bulk.'+str(bulkStart).rjust(7,'0')+'.vlsv'
print(bulkname)

f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    # Get size of simulation (in cells)
xsize = f.read_parameter('xcells_ini')
#zsize = f.read_parameter('zcells_ini')
ysize = f.read_parameter('ycells_ini')

    # Spatial size
xmax = f.read_parameter("xmax")
xmin = f.read_parameter("xmin")
dx   = (xmax-xmin)/xsize

#zmax = f.read_parameter("zmax")
#zmin = f.read_parameter("zmin")
#dz   = (zmax-zmin)/zsize

ymax = f.read_parameter("ymax")
ymin = f.read_parameter("ymin")
dy   = (ymax-ymin)/ysize

    # Box coordinates # m
RE = 6371e3 # m
x1 = int((3.0 *RE - xmin) / dx)
x2 = int((6.0 *RE - xmin) / dx)
#z1 = int((11.0 *RE - zmin) / dz)
#z2 = int((14.0 *RE - zmin) / dz)
y1 = int((11.0 *RE - ymin) / dy)
y2 = int((14.0 *RE - ymin) / dy)


# One point in the middle
xmid = int((x2 + x1) / 2)
ymid = int((y2 + y1) / 2)

# Initialise field
Bmid = np.zeros([3,bulkTot])

print('Start processing files')
for i in range(0,bulkTot):
    bulkname  = 'bulk.'+str(i+bulkStart).rjust(7,'0')+'.vlsv'
    print(bulkname)

    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    # Read cellids. Cellids are randomly sorted because of multiprocessing
    cellids = f.read_variable('CellID')

    # Read field
    B = f.read_variable('B')
    # Sort the field to be in the order of cellids (1 to N). Reshape it to be in the 2D simulation frame.
    B         = B[cellids.argsort()].reshape([ysize,xsize,3]) # Starts from bottom left (y<0,x<0) to the right
    Bmid[:,i] = B[ymid,xmid,:]
print('All files processed')


# Time avg of B per component & norm
B_avg  = np.mean(Bmid,axis=1) 
B_norm = np.linalg.norm(B_avg)

# Angle
ez = [0,0,1]
angle = abs(np.arccos(np.dot(B_avg,ez)/B_norm) * 180.0 / np.pi -180.0) # Degrees

np.save(path_save+'angle.npy',angle)
print('Angle = '+str(angle)+', saved.')
