import pytools as pt
import numpy as np
import scipy
import scipy.ndimage
import sys

path_bulk = '/proj/vlasov/2D/ABA/bulk/'
path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/ABA/Data/FFT/'

angle_rot = - np.load(path_save+'angle.npy')
print angle_rot

if len(sys.argv)==3: # Starting and end bulk given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting bulk given, generate one bulk
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)
  
    # Open file
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
    #z1 = int((15.0 *RE - zmin) / dz)
    #z2 = int((18.0 *RE - zmin) / dz)   
    y1 = int((11.0 *RE - ymin) / dy)
    y2 = int((14.0 *RE - ymin) / dy)   

    # Read cellids. Cellids are randomly sorted because of multiprocessing
    cellids = f.read_variable('CellID')
    # Read field
    B  = f.read_variable('B')
    E  = f.read_variable('E')
    # Get right component
    Bx = B[:,0]
    By = B[:,1]
    Bz = B[:,2]
    Ex = E[:,0]
    Ey = E[:,1]
    Ez = E[:,2]
    # Sort the field to be in the order of cellids (1 to N). Reshape it to be in the 2D simulation frame.
    Bx = Bx[cellids.argsort()].reshape([ysize,xsize]) #Starts from bottom left (y<0,x<0) to the right
    Bx = Bx[y1:y2,x1:x2]                             #Selected box. Comment if whole simulation needed.
    By = By[cellids.argsort()].reshape([ysize,xsize]) 
    By = By[y1:y2,x1:x2]
    Bz = Bz[cellids.argsort()].reshape([ysize,xsize]) 
    Bz = Bz[y1:y2,x1:x2]
    Ex = Ex[cellids.argsort()].reshape([ysize,xsize]) 
    Ex = Ex[y1:y2,x1:x2]                             
    Ey = Ey[cellids.argsort()].reshape([ysize,xsize])     
    Ey = Ey[y1:y2,x1:x2]
    Ez = Ez[cellids.argsort()].reshape([ysize,xsize]) 
    Ez = Ez[y1:y2,x1:x2]
    # Rotate for FFT
    BxR = scipy.ndimage.rotate(Bx,angle_rot,reshape=False)
    ByR = scipy.ndimage.rotate(By,angle_rot,reshape=False)
    BzR = scipy.ndimage.rotate(Bz,angle_rot,reshape=False)
    ExR = scipy.ndimage.rotate(Ex,angle_rot,reshape=False)
    EyR = scipy.ndimage.rotate(Ey,angle_rot,reshape=False)
    EzR = scipy.ndimage.rotate(Ez,angle_rot,reshape=False)
    #Project to 1D. Will sum over all waves in the box. Might add some features if waves are not in same direction as cut.
    #F1D = np.zeros(Field.shape[1]) #Sum over y-axis. Switch to Field.shape[0] to sum over x-axis.
    #F1D = np.sum(Field, axis = 0)  #Sum over y-axis. Switch to axis = 1 to sum over x-axis.

    # K_PERP
        # Slice to 1D. Will take values over a line in the middle of the box. Might see waves moving faster than they are.
    Bx1D_kperp = BxR[BxR.shape[0]/2,:] 
    By1D_kperp = ByR[ByR.shape[0]/2,:]
    Bz1D_kperp = BzR[BzR.shape[0]/2,:]
    Ex1D_kperp = ExR[ExR.shape[0]/2,:]
    Ey1D_kperp = EyR[EyR.shape[0]/2,:]
    Ez1D_kperp = EzR[EzR.shape[0]/2,:]
        # Hanning windowing to force the edge of box to be periodic. Uncomment if boundary conditions arent periodic.
    Bx1D_kperp *= np.hanning(len(Bx1D_kperp)) #Field.shape[1] if along y-axis, len(F1D) if domain is a line
    By1D_kperp *= np.hanning(len(By1D_kperp))
    Bz1D_kperp *= np.hanning(len(Bz1D_kperp))
    Ex1D_kperp *= np.hanning(len(Ex1D_kperp))
    Ey1D_kperp *= np.hanning(len(Ey1D_kperp))
    Ez1D_kperp *= np.hanning(len(Ez1D_kperp))
        # Save
    np.save(path_save+'Bx1D_kperp_'+str(j)+'.npy',Bx1D_kperp)
    np.save(path_save+'By1D_kperp_'+str(j)+'.npy',By1D_kperp)
    np.save(path_save+'Bz1D_kperp_'+str(j)+'.npy',Bz1D_kperp)
    np.save(path_save+'Ex1D_kperp_'+str(j)+'.npy',Ex1D_kperp)
    np.save(path_save+'Ey1D_kperp_'+str(j)+'.npy',Ey1D_kperp)
    np.save(path_save+'Ez1D_kperp_'+str(j)+'.npy',Ez1D_kperp)

    # K_PARA
        # Slice to 1D. Will take values over a line in the middle of the box. Might see waves moving faster than they are.
    Bx1D_kpara = BxR[:,BxR.shape[1]/2]
    By1D_kpara = ByR[:,ByR.shape[1]/2]
    Bz1D_kpara = BzR[:,BzR.shape[1]/2]
    Ex1D_kpara = ExR[:,ExR.shape[1]/2]
    Ey1D_kpara = EyR[:,EyR.shape[1]/2]
    Ez1D_kpara = EzR[:,EzR.shape[1]/2]
        # Hanning windowing to force the edge of box to be periodic. Uncomment if boundary conditions arent periodic.
    Bx1D_kpara *= np.hanning(len(Bx1D_kpara)) #Field.shape[1] if along y-axis, len(F1D) if domain is a line
    By1D_kpara *= np.hanning(len(By1D_kpara))
    Bz1D_kpara *= np.hanning(len(Bz1D_kpara))
    Ex1D_kpara *= np.hanning(len(Ex1D_kpara))
    Ey1D_kpara *= np.hanning(len(Ey1D_kpara))
    Ez1D_kpara *= np.hanning(len(Ez1D_kpara))
        # Save
    np.save(path_save+'Bx1D_kpara_'+str(j)+'.npy',Bx1D_kpara)
    np.save(path_save+'By1D_kpara_'+str(j)+'.npy',By1D_kpara)
    np.save(path_save+'Bz1D_kpara_'+str(j)+'.npy',Bz1D_kpara)
    np.save(path_save+'Ex1D_kpara_'+str(j)+'.npy',Ex1D_kpara)
    np.save(path_save+'Ey1D_kpara_'+str(j)+'.npy',Ey1D_kpara)
    np.save(path_save+'Ez1D_kpara_'+str(j)+'.npy',Ez1D_kpara)
    
