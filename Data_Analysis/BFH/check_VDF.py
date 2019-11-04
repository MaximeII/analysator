# Get data and plot VDF for a cell in the middle of square where FFT is taken

import pytools as pt
import numpy as np
import pylab
import sys

path_bulk = '/scratch/project_2000203/2D/BFH/reverted_ionosphere_field_boundary/'
path_save = '/users/dubartma/analysator/Data_Analysis/BFH/Data/'

RE = 6371e3 # m
coord_mid = [(3.0 + 6.0)/2 *RE,0.0,(15.0 + 18.0)/2.0 *RE]

# Qpara Nose: 0.0 3.0 | -16.5 -13.5
# Qpara Tail: -31.5 -28.5 | -46.5 -43.5
# Qperp Nose: 3.0 6.0 | 15.0 18.0
# Qperp Tail: -31.5 -28.5 | 28.5 31.5
# BCH Nose: 0.0 3.0 | 20.0 23.0
# BCH Tail: -31.5 -28.5 | 38.5 41.5

if len(sys.argv)==3: # Starting and end bulk given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting bulk given, generate one bulk
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)
    #distribname = 'distributions.'+str(int(j/10)).rjust(7,'0')+'.vlsv'
    #print(distribname)

    #Open file
    #f = pt.vlsvfile.VlsvReader(path_VDF+distribname)
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    CellID = int(f.get_cellid(coord_mid))
    print(CellID)

    cell_candidates = f.read(mesh = "SpatialGrid", tag = "CELLSWITHBLOCKS")
    # Read in the coordinates of the cells:
    cell_candidate_coordinates = [f.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]

    # Read in the cell's coordinates:
    pick_cell_coordinates = f.get_cell_coordinates(CellID)
    if len(cell_candidates) == 0:
        print("No velocity distribution data found in this file!")
        sys.exit()

    # Find the nearest:
    from operator import itemgetter
    norms = np.sum((cell_candidate_coordinates - pick_cell_coordinates)**2, axis=-1)**(1./2)
    norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))

    # Get the cell id:
    cid = cell_candidates[i]
    print("PLOTTED CELL ID: " + str(cid))

    np.save(path_save+'CellID_SC.npy',cid)

    # check if velocity space exists in this cell
    if f.check_variable('fSaved'): #restart files will not have this value                                                                      
        if f.read_variable('fSaved',cid) != 1.0:
            print('No vdf in this cell')

    print("CHECKED")
 
    #pt.plot.plot_vdf(filename=path_bulk+bulkname,cellids=cid,wmark=None,yz=True)
    #pylab.plt.savefig(path_fig+'yz/VDF_'+str(j).rjust(7,'0')+'.png',dpi=400)
    #print(path_fig+'xy/VDF_'+str(j).rjust(7,'0')+'.png')

