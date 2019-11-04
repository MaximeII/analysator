#Does 2D FFT of provided data.

import pytools as pt
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import imageio
import scipy
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

path_HYDROS = '/homeappl/home/dubartma/HYDROS/HYDROS/output/'
run       = 'BCQ'
path_bulk = '/proj/vlasov/2D/'+run+'/bulk/'
bulkStart = 1600
bulkEnd   = 2400
bulkTot   = bulkEnd - bulkStart + 1
#path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/'+run+'/Data/FFT/'
#path_fig  = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/'+run+'/Fig/FFT/'
path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/Data/Qperp/Nose/FFT/'
path_fig  = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/Fig/Qperp/Nose/FFT/'
k_str = ['kpara','kperp']

len_y = np.zeros(len(k_str))
for i in range(0,len(k_str)):
    len_y[i] = len(np.load(path_save+'Bx1D_'+k_str[i]+'_'+str(bulkStart)+'.npy'))

try:
    Bx_kpara = np.load(path_save+'Bx_'+k_str[0]+'.npy')
    By_kpara = np.load(path_save+'By_'+k_str[0]+'.npy')
    Bz_kpara = np.load(path_save+'Bz_'+k_str[0]+'.npy')
    Ex_kpara = np.load(path_save+'Ex_'+k_str[0]+'.npy')
    Ey_kpara = np.load(path_save+'Ey_'+k_str[0]+'.npy')
    Ez_kpara = np.load(path_save+'Ez_'+k_str[0]+'.npy')
    Bx_kperp = np.load(path_save+'Bx_'+k_str[1]+'.npy')
    By_kperp = np.load(path_save+'By_'+k_str[1]+'.npy')
    Bz_kperp = np.load(path_save+'Bz_'+k_str[1]+'.npy')
    Ex_kperp = np.load(path_save+'Ex_'+k_str[1]+'.npy')
    Ey_kperp = np.load(path_save+'Ey_'+k_str[1]+'.npy')
    Ez_kperp = np.load(path_save+'Ez_'+k_str[1]+'.npy')
except IOError:

    Bx_kpara = np.zeros([len_y[0],bulkTot])
    By_kpara = np.zeros([len_y[0],bulkTot])
    Bz_kpara = np.zeros([len_y[0],bulkTot])
    Ex_kpara = np.zeros([len_y[0],bulkTot])
    Ey_kpara = np.zeros([len_y[0],bulkTot])
    Ez_kpara = np.zeros([len_y[0],bulkTot])
    Bx_kperp = np.zeros([len_y[1],bulkTot])
    By_kperp = np.zeros([len_y[1],bulkTot])
    Bz_kperp = np.zeros([len_y[1],bulkTot])
    Ex_kperp = np.zeros([len_y[1],bulkTot])
    Ey_kperp = np.zeros([len_y[1],bulkTot])
    Ez_kperp = np.zeros([len_y[1],bulkTot])
    for i in range(0,bulkTot):
        Bx_kpara[:,i] = np.load(path_save+'Bx1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        By_kpara[:,i] = np.load(path_save+'By1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Bz_kpara[:,i] = np.load(path_save+'Bz1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ex_kpara[:,i] = np.load(path_save+'Ex1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ey_kpara[:,i] = np.load(path_save+'Ey1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ez_kpara[:,i] = np.load(path_save+'Ez1D_'+k_str[0]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Bx_kperp[:,i] = np.load(path_save+'Bx1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        By_kperp[:,i] = np.load(path_save+'By1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Bz_kperp[:,i] = np.load(path_save+'Bz1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ex_kperp[:,i] = np.load(path_save+'Ex1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ey_kperp[:,i] = np.load(path_save+'Ey1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]
        Ez_kperp[:,i] = np.load(path_save+'Ez1D_'+k_str[1]+'_'+str(i+bulkStart)+'.npy') * np.hanning(bulkTot)[i]

    np.save(path_save+'Bx_'+k_str[0]+'.npy',Bx_kpara)
    np.save(path_save+'By_'+k_str[0]+'.npy',By_kpara)
    np.save(path_save+'Bz_'+k_str[0]+'.npy',Bz_kpara)
    np.save(path_save+'Ex_'+k_str[0]+'.npy',Ex_kpara)
    np.save(path_save+'Ey_'+k_str[0]+'.npy',Ey_kpara)
    np.save(path_save+'Ez_'+k_str[0]+'.npy',Ez_kpara)
    np.save(path_save+'Bx_'+k_str[1]+'.npy',Bx_kperp)
    np.save(path_save+'By_'+k_str[1]+'.npy',By_kperp)
    np.save(path_save+'Bz_'+k_str[1]+'.npy',Bz_kperp)
    np.save(path_save+'Ex_'+k_str[1]+'.npy',Ex_kperp)
    np.save(path_save+'Ey_'+k_str[1]+'.npy',Ey_kperp)
    np.save(path_save+'Ez_'+k_str[1]+'.npy',Ez_kperp)

print("Arrays' building done")

#2D Fourier Transform
k_omega_Bx_kpara = np.fft.fftshift(np.fft.fft2(Bx_kpara.T)) #Shift to be centered around 0. Transposed to have complex part.
k_omega_By_kpara = np.fft.fftshift(np.fft.fft2(By_kpara.T))
k_omega_Bz_kpara = np.fft.fftshift(np.fft.fft2(Bz_kpara.T))
k_omega_Ex_kpara = np.fft.fftshift(np.fft.fft2(Ex_kpara.T))
k_omega_Ey_kpara = np.fft.fftshift(np.fft.fft2(Ey_kpara.T))
k_omega_Ez_kpara = np.fft.fftshift(np.fft.fft2(Ez_kpara.T))
k_omega_Bx_kperp = np.fft.fftshift(np.fft.fft2(Bx_kperp.T)) #Shift to be centered around 0. Transposed to have complex part.
k_omega_By_kperp = np.fft.fftshift(np.fft.fft2(By_kperp.T))
k_omega_Bz_kperp = np.fft.fftshift(np.fft.fft2(Bz_kperp.T))
k_omega_Ex_kperp = np.fft.fftshift(np.fft.fft2(Ex_kperp.T))
k_omega_Ey_kperp = np.fft.fftshift(np.fft.fft2(Ey_kperp.T))
k_omega_Ez_kperp = np.fft.fftshift(np.fft.fft2(Ez_kperp.T))
print '2D FFT Done'

#Get size of simulation (in cells)
f       = pt.vlsvfile.VlsvReader(path_bulk+'bulk.'+str((bulkStart+bulkEnd)/2).rjust(7,'0')+'.vlsv')
xsize   = f.read_parameter('xcells_ini')
zsize   = f.read_parameter('zcells_ini')
#ysize = f.read_parameter('ycells_ini')
cellids = f.read_variable("CellID")

#Spatial size
xmax = f.read_parameter("xmax")
xmin = f.read_parameter("xmin")
dx   = (xmax-xmin)/xsize

zmax = f.read_parameter("zmax")
zmin = f.read_parameter("zmin")
dz   = (zmax-zmin)/zsize

#ymax = f.read_parameter("ymax")
#ymin = f.read_parameter("ymin")
#dy   = (ymax-ymin)/ysize


#Box coordinates # m
RE = 6371e3 # m
x1 = int((3.0 *RE - xmin) / dx)
x2 = int((6.0 *RE - xmin) / dx)
z1 = int((15.0 *RE - zmin) / dz)
z2 = int((18.0 *RE - zmin) / dz)
#y1 = int((15.0 *RE - ymin) / dy)
#y2 = int((18.0 *RE - ymin) / dy)

# Qpara Nose: 0.0 3.0 | -16.5 -13.5
# Qpara Tail: -31.5 -28.5 | -46.5 -43.5
# Qperp Nose: 3.0 6.0 | 15.0 18.0
# Qperp Tail: -31.5 -28.5 | 28.5 31.5
# BCH Nose: 0.0 3.0 | 20.0 23.0
# BCH Tail: -31.5 -28.5 | 38.5 41.5

#Simulation time step
dt = 0.5

#Speed
v = f.read_variable("v")
v = v[cellids.argsort()].reshape([zsize,xsize,3])
v = v[z1:z2,x1:x2,:]
#vnorm = np.sqrt(v[:,:,0]**2 + v[:,:,1]**2 + v[:,:,2]**2)
#vnorm = np.mean(vnorm)
v = [np.mean(v[:,:,0]),np.mean(v[:,:,1]),np.mean(v[:,:,2])]
vnorm = np.linalg.norm(v)

#Bfield
B     = f.read_variable('B')
B     = B[cellids.argsort()].reshape([zsize,xsize,3])
B     = B[z1:z2,x1:x2,:]
B     = [np.mean(B[:,:,0]),np.mean(B[:,:,1]),np.mean(B[:,:,2])]
Bnorm = np.linalg.norm(B)

vpara = - np.dot(B,v)/Bnorm
vperp = - np.linalg.norm(np.cross(B,np.cross(B,v))/Bnorm**2)

#Density
rho = f.read_variable('rho')
rho = rho[cellids.argsort()].reshape([zsize,xsize])
rho = rho[z1:z2,x1:x2]
rho = np.mean(rho)

#Pressure
P     = f.read_variable('Pressure')
P     = P[cellids.argsort()].reshape([zsize,xsize])
P     = P[z1:z2,x1:x2]
Pnorm = np.mean(P)

#Constants (SI units)
q        = 1.60217e-19
mu0      = 1.256637061e-6
c        = 299792458
mp       = 1.6726219e-27
kB       = 1.38064852e-23
me       = 9.1093826e-31
epsilon0 = 8.85418781762e-12
gamma    = 5.0/3.0
d_i      = np.sqrt(mp*epsilon0*c*c/(rho*q*q))
w_ci     = q*Bnorm/mp
vA       = Bnorm / np.sqrt(mu0 * rho * mp)
vS       = np.sqrt(gamma * Pnorm / (rho * mp))


#Dispersion relations
def CFL(k): #CFL condition
    return dx/0.044304 * (k / d_i)   #Constante may change with simulation

def AWaves_para(k): #Alfven waves
    return abs(Bnorm / np.sqrt(mu0 * rho * mp) * (abs(-k)/d_i) + vpara*(-k/d_i)) # - flow speed

def AWaves_perp(k): #Alfven waves
    return abs( Bnorm / np.sqrt(mu0 * rho * mp) * (abs(-k)/d_i) - vperp*(-k/d_i))

def FMW_kpara(k): #Fast Magnetosonic Waves
    return abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 + np.sqrt( (vS**2 + vA**2)**2 - 4* vA**2 * vS**2 ))) * (abs(-k)/d_i) + vpara*(-k/d_i))

def SMW_kpara(k): #Slow Magnetosonic Waves
    return abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 - np.sqrt( (vS**2 + vA**2)**2 - 4* vA**2 * vS**2 ))) * (abs(-k)/d_i) + vpara*(-k/d_i))

def Bulk_v_para(k): #Bulk Velocity
    return  abs(vpara*(-k/d_i))

def FMW_kperp(k): #Fast Magnetosonic Waves
    return  abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 + np.sqrt( (vS**2 + vA**2)**2))) * (abs(-k)/d_i) - vperp*(-k/d_i))

def SMW_kperp(k): #Slow Magnetosonic Waves
    return  abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 - np.sqrt( (vS**2 + vA**2)**2))) * (abs(-k)/d_i) - vperp*(-k/d_i))

def Bulk_v_perp(k): #Bulk Velocity
    return  abs(vperp*(-k/d_i))

#Plot Fourier space
kmin     =- np.pi / dx * d_i #Change to dy/dz if needed
kmax     = np.pi / dx * d_i
omegamin = 0
omegamax = 2.8 #np.pi / dt / w_ci

kmin_plot = - np.pi / 300000.0 * d_i
kmax_plot = np.pi / 300000.0 * d_i

a      = np.linspace(kmin, kmax, 2000)
cfl    = map(abs,map(CFL,a))/w_ci
awaves_para = map(AWaves_para,a)/w_ci
awaves_perp = map(AWaves_perp,a)/w_ci
fmw_kperp = map(FMW_kperp,a)/w_ci
smw_kperp = map(SMW_kperp,a)/w_ci
fmw_kpara = map(FMW_kpara,a)/w_ci
smw_kpara = map(SMW_kpara,a)/w_ci
bv_para   = map(Bulk_v_para,a)/w_ci
bv_perp   = map(Bulk_v_perp,a)/w_ci

Bmax = 1e-4
Bmin = 1e-13
Emax = 1e1
Emin = 1e-7

fig, (ax_para, ax_perp) = plt.subplots(1,2,sharey=True)
labelsize = 14
fontsize  = 15

#K_PARA
image = ax_para.imshow(abs(k_omega_Ey_kpara[(k_omega_Ey_kpara.shape[0])/2:,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')

ax_para.locator_params(axis='x', nbins=5)      #Set number of ticks
ax_para.tick_params(axis='both',labelsize=labelsize)   #Set fontsize of tick labels
ax_para.xaxis.set_major_locator(MultipleLocator(0.5))

ax_para.set_xlim(kmin_plot, kmax_plot)
ax_para.set_ylim(omegamin, omegamax)
ax_para.set_xlabel("$k_\\parallel * d_i$",fontsize=fontsize)
ax_para.set_ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize=fontsize)

ax_para.plot(a, cfl, linewidth=2, color='black',label='CFL condition')
ax_para.plot(a,awaves_para,linewidth=2,color='blue',label='Alfven Waves')
ax_para.plot(a,fmw_kpara,linewidth=2,color='red',linestyle='--',label='Fast MW')
ax_para.plot(a,smw_kpara,linewidth=2,color='red',linestyle=':',label='Slow MW')
ax_para.plot(a,bv_para,linewidth=2,color='black',linestyle='--',label='V')

ax_para.text(kmax_plot - 0.3,2.7,'(a)',fontsize=12)

ax_para.legend(bbox_to_anchor = [-0.03, 1.04, 2.11, 0.05], loc='upper right',ncol=5, mode='expand', fontsize=11)

#K_PERP
image = ax_perp.imshow(abs(k_omega_Ey_kperp[(k_omega_Ey_kperp.shape[0])/2:,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')

ax_perp.locator_params(axis='x', nbins=5) #Set number of ticks
ax_perp.tick_params(axis='x',labelsize=labelsize) #Set fontsize of tick labels
ax_perp.xaxis.set_major_locator(MultipleLocator(0.5))

ax_perp.set_xlim(kmin_plot, kmax_plot)
ax_perp.set_ylim(omegamin, omegamax)
ax_perp.set_xlabel("$k_\\bot * d_i$",fontsize=fontsize)

ax_perp.plot(a, cfl, linewidth=2, color='black',label='CFL condition')
ax_perp.plot(a,awaves_perp,linewidth=2,color='blue',label='Alfven Speed')
ax_perp.plot(a,fmw_kperp,linewidth=2,color='red',linestyle='--',label='Fast MW')
ax_perp.plot(a,smw_kperp,linewidth=2,color='red',linestyle=':',label='Slow MW')
ax_perp.plot(a,bv_perp,linewidth=2,color='black',linestyle='--',label='V')

ax_perp.text(kmax_plot - 0.3,2.7,'(b)',fontsize=12)

fig.suptitle('$\\Delta r = 300$ km',fontsize=16)

plt.tight_layout(rect=[-0.01, 0.0, 1, 0.92])

#Colorbar
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.86, 0.116, 0.03, 0.775])
fig.colorbar(image, cax=cbar_ax)
cbar_ax.set_ylabel('$E_y^2 (V/m)^2$',fontsize=fontsize)


title = 'FFT_'+run 
#plt.savefig(path_fig+title+'.eps',dpi=400)
plt.savefig(path_fig+title+'.pdf',dpi=400)
print(path_fig+title+'.pdf')


