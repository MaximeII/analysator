import pytools as pt
import sys, os, socket
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import matplotlib
import matplotlib.patches as patches

# Avoids opening a figure window
if str(matplotlib.get_backend()) is not 'Agg':
    plt.switch_backend('Agg') 

run = 'BCQ'
fileLocation='/proj/vlasov/2D/'+run+'/bulk/'
fluxLocation='/proj/vlasov/2D/'+run+'/flux/'
outputLocation=outputdir='/homeappl/home/dubartma/appl_taito/analysator/Whistler/'+run+'/'
#path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/ABA/Data/'
path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/'+run+'/Data/'

# Frame extent for this job given as command-line arguments
if len(sys.argv)==3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    # CellIDs of interest
    cid1 = int(np.load(path_save+'CellID_SC.npy'))
    Re = 6371e3 # m

    # Figure layout
    mh = 0.05
    mv = 0.12
    ecb = 0.015
    lcb = 0.02
    ev = 0.1
    h1 = 1-1.5*mv
    h2 = (h1-ev)/2.
    l1 = 0.32
    l2 = l1*h2/h1
    eh = (1.-l1-2*l2-2.5*ecb-2*lcb-2*mh)/3.
    fontscale = 2.5
    figratio = h2/l2

    fig = plt.figure(figsize=(5.*figratio,5.))
    # Colormap and its colorbar
    ax_col = fig.add_axes([mh,mv,l1,h1])
    cax_col = fig.add_axes([mh+l1+ecb,mv,lcb,h1])
    # VDFs and their colorbar
    deltaX = mh+l1+ecb+lcb+2*eh
    ax_VDF1_par  = fig.add_axes([deltaX,mv+h2+ev,l2,h2])
    ax_VDF1_per  = fig.add_axes([deltaX+l2+eh,mv+h2+ev,l2,h2])
    ax_VDF1_per2 = fig.add_axes([deltaX,mv,l2,h2])
   # ax_VDF2_per = fig.add_axes([deltaX+l2+eh,mv,l2,h2])
    cax_VDF = fig.add_axes([deltaX+2*l2+eh+1.5*ecb,mv,lcb,h1-0.01])


    # Colormap with J_diamagnetic
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var='B',
                          vmin=1.0e-8,vmax=3.0e-8,
                          boxre=[0.0,9.0,12.0,21.0],
                          run="BCQ",
                          colormap='hot_desaturated',
                          scale=1.3,
                          axes=ax_col,cbaxes=cax_col,
                          outputdir=outputLocation,lin=1)

    rect = patches.Rectangle((3.0,15.0),3.0,3.0,linewidth=2,edgecolor='red',facecolor='none')
    ax_col.add_patch(rect)

    # Addition of location of cells with VDFs
    vlsvReader=pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    xCid,yCid,zCid = vlsvReader.get_cell_coordinates(cid1)
    VDFcoord_1 = [xCid/Re,yCid/Re,zCid/Re]

    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='black',marker='o',s=20)
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='white',marker='o',s=2)

    # VDFs of first cell
    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bpara1=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_par,cbaxes=cax_VDF,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=1.3,
                     colormap='nipy_spectral')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bperp=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=1.3,
                     colormap='nipy_spectral')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bpara=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per2,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=1.3,
                     colormap='nipy_spectral')

   # # VDFs of second cell
   # pt.plot.plot_vdf(filename=fileLocation+bulkname,
   #                  cellids=[cid2],
   #                  bpara1=1,
   #                  box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
   #                  fmin=1e-15, fmax=1e-11,
   #                  slicethick=1,
   #                  axes=ax_VDF2_par,nocb=1,
   #                  cbulk=1,
   #                  axisunit=6,
   #                  title='',
   #                  scale=1.15,
   #                  colormap='nipy_spectral')

   # pt.plot.plot_vdf(filename=fileLocation+bulkname,
   #                  cellids=[cid2],
   #                  bperp=1,
   #                  box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
   #                  fmin=1e-15, fmax=1e-11,
   #                  slicethick=1,
   #                  axes=ax_VDF2_per,nocb=1,
   #                  cbulk=1,
   #                  axisunit=6,
   #                  title='',
   #                  scale=1.15,
   #                  colormap='nipy_spectral')
    fig.suptitle('(a) $\\Delta r = 300$ km',fontsize=15)
    figname=run+"_VDFs_"+str(j).rjust(7,'0')
    extension = '.png'
    plt.savefig(outputLocation+figname+extension,dpi=800)
    print(outputLocation+figname+extension)
    #plt.savefig(outputLocation+figname+'.eps')
