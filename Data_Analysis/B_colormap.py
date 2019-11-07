import pytools as pt
import numpy as np 
import matplotlib.pyplot as plt
import sys
from matplotlib import rcParams

outputLocation=outputdir='/users/dubartma/analysator/Data_Analysis/Fig/'
hspace = 0.02
scale  = 1.4

if len(sys.argv)==3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    fig = plt.figure(figsize = (10.,5.))
   
    run = 'BCQ'
    fileLocation='/scratch/project_2000203/2D/'+run+'/'

    title = '(a) t = '+str(j/2.0)+' s, $\\Delta r = 300$ km'

    x1 = hspace
    y1 = 0.1
    h1 = 0.76
    l1 = 0.32

    ax_col1 = fig.add_axes([x1,y1,l1,h1])

    # Colormap of Temperature anisotropy for BCQ
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var='B',
                          vmin=1.0e-8,vmax=3.0e-8,
                          boxre=[0.0,17.0,-30.0,30.0],
                          run="BCQ",
                          #title=title,
                          colormap='hot_desaturated',lin=1,
                          scale=scale,
                          axes=ax_col1,nocb=1,
                          outputdir=outputLocation)
  
    title1 = ax_col1.set_title(title,fontsize=15)
    title1.set_position([0.5,1.09])    


    run = 'BFH'
    fileLocation = '/scratch/project_2000203/2D/'+run+'/reverted_ionosphere_field_boundary/'
    bulk = 1000
    bulknameBFH = 'bulk.0001000.vlsv'

    title = '(b) t = '+str(bulk/2.0)+' s, $\\Delta r = 600$ km' 

    x2 = 2.0*hspace + 1.0*l1
    y2 = y1
    h2 = h1
    l2 = l1

    ax_col2  = fig.add_axes([x2,y2,l2,h2])

    # Colormap of Temperature anisotropy for BFH
    
    pt.plot.plot_colormap(filename=fileLocation+bulknameBFH,var='fg_b',
                          vmin=1.0e-8,vmax=3.0e-8,
                          boxre=[0.0,17.0,-30.0,30.0],
                          run="BFG",
                          #title=title,
                          colormap='hot_desaturated',lin=1,
                          scale=scale,
                          axes=ax_col2,nocb=1,
                          outputdir=outputLocation)


    title2 = ax_col2.set_title(title,fontsize=15)
    title2.set_position([0.5,1.09])


    run = 'BCG'
    fileLocation='/scratch/project_2000203/2D/'+run+'/'

    title = '(c) t = '+str(j/2.0)+' s, $\\Delta r = 900$ km'

    x3 = 3.0*hspace + 2.0*l1
    y3 = y1
    h3 = h1
    l3 = l1

    xcb = x3+l3+0.015
    ycb = y3
    lcb = 0.05
    hcb = h1

    ax_col3  = fig.add_axes([x3,y3,l3,h3])

    # Colormap of Temperature anisotropy for BCG
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var='B',
                          vmin=1.0e-8,vmax=3.0e-8,
                          boxre=[0.0,17.0,-30.0,30.0],
                          run="BCG",
                          #title=title,
                          colormap='hot_desaturated',lin=1,
                          scale=scale,
                          axes=ax_col3,
                          outputdir=outputLocation) 

    title3 = ax_col3.set_title(title,fontsize=15,pad=30)
    #title3.set_position([0.5,1.09])

    figname = 'B_colormap'
    plt.savefig(outputLocation+figname+'.png',dpi=800)
    print(outputLocation+figname+'.png')
