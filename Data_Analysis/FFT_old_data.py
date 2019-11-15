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
run       = 'BCG'
path_bulk = '/scratch/project_2000203/2D/'+run+'/'
path_save = '/users/dubartma/analysator/Data_Analysis/Data/'
outputdir = '/users/dubartma/analysator/Data_Analysis/Fig/'

k_omega_Ey_kpara = np.load(path_save+'k_omega_Ey_kpara_'+run+'.npy')
k_omega_Ey_kperp = np.load(path_save+'k_omega_Ey_kperp_'+run+'.npy')

#Get size of simulation (in cells)
f       = pt.vlsvfile.VlsvReader(path_bulk+'bulk.0002000.vlsv')
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

#Speed
v = f.read_variable("v")
v = v[cellids.argsort()].reshape([zsize,xsize,3])
v = v[z1:z2,x1:x2,:]
#vnorm = np.sqrt(v[:,:,0]**2 + v[:,:,1]**2 + v[:,:,2]**2)
#vnorm = np.mean(vnorm)
v = [np.mean(v[:,:,0]),np.mean(v[:,:,1]),np.mean(v[:,:,2])]
vnorm = np.linalg.norm(v)
print('Vnorm = '+str(vnorm)+' m/s')

#Bfield
B     = f.read_variable('B')
B     = B[cellids.argsort()].reshape([zsize,xsize,3])
B     = B[z1:z2,x1:x2,:]
B     = [np.mean(B[:,:,0]),np.mean(B[:,:,1]),np.mean(B[:,:,2])]
Bnorm = np.linalg.norm(B)
print('Bnorm = '+str(Bnorm)+' T')

vpara = - np.dot(B,v)/Bnorm
vperp = - np.linalg.norm(np.cross(B,np.cross(B,v))/Bnorm**2)
print('Vpara = '+str(vpara)+' m/s')
print('Vperp = '+str(vperp)+' m/s')

#Density
rho = f.read_variable('rho')
rho = rho[cellids.argsort()].reshape([zsize,xsize])
rho = rho[z1:z2,x1:x2]
rho = np.mean(rho)
print('n = '+str(rho)+' m^-3')

#Pressure
P     = f.read_variable('Pressure')
P     = P[cellids.argsort()].reshape([zsize,xsize])
P     = P[z1:z2,x1:x2]
Pnorm = np.mean(P)
print('Pnorm = '+str(Pnorm)+' Pa')

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

print('d_i = '+str(d_i)+' m')
print('w_ci = '+str(w_ci)+' s^-1')
print('vA = '+str(vA)+' m/s')
print('vS = '+str(vS)+' m/s')

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
omegamax = 2.8 #np.pi / dt / w_c

kmin_plot = - np.pi / 300000.0 * d_i
kmax_plot = np.pi / 300000.0 * d_i

print(kmin_plot,kmax_plot)

a           = np.linspace(kmin, kmax, 2000)
cfl         = list(map(abs,map(CFL,a)))/w_ci
awaves_para = list(map(AWaves_para,a))/w_ci
awaves_perp = list(map(AWaves_perp,a))/w_ci
fmw_kperp   = list(map(FMW_kperp,a))/w_ci
smw_kperp   = list(map(SMW_kperp,a))/w_ci
fmw_kpara   = list(map(FMW_kpara,a))/w_ci
smw_kpara   = list(map(SMW_kpara,a))/w_ci
bv_para     = list(map(Bulk_v_para,a))/w_ci
bv_perp     = list(map(Bulk_v_perp,a))/w_ci

Emax = 1e1
Emin = 1e-7

fig, (ax_para, ax_perp) = plt.subplots(1,2,sharey=True)
labelsize = 11
fontsize  = 15
linewidth = 1

cmap = 'jet'
#K_PARA
image = ax_para.imshow(abs(k_omega_Ey_kpara[int((k_omega_Ey_kpara.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto', cmap = plt.get_cmap(cmap))

ax_para.locator_params(axis='x', nbins=5)      #Set number of ticks
ax_para.tick_params(axis='both',labelsize=labelsize)   #Set fontsize of tick labels
ax_para.xaxis.set_major_locator(MultipleLocator(0.5))

ax_para.set_xlim(kmin_plot, kmax_plot)
ax_para.set_ylim(omegamin, omegamax)
ax_para.set_xlabel("$k_\\parallel * d_i$",fontsize=fontsize)
ax_para.set_ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize=fontsize)

ax_para.plot(a, cfl, linewidth=linewidth, color='black',label='CFL condition')
ax_para.plot(a,awaves_para,linewidth=linewidth,color='blue',label='Alfven Waves')
ax_para.plot(a,fmw_kpara,linewidth=linewidth,color='red',linestyle='--',label='Fast MW')
ax_para.plot(a,smw_kpara,linewidth=linewidth,color='red',linestyle=':',label='Slow MW')
ax_para.plot(a,bv_para,linewidth=linewidth,color='black',linestyle='--',label='V')

ax_para.text(kmax_plot - 0.3,2.7,'(e)',fontsize=12)

#ax_para.legend(bbox_to_anchor = [-0.03, 1.04, 2.105, 0.05], loc='upper right',ncol=5, mode='expand', fontsize=8.5)

#K_PERP
image = ax_perp.imshow(abs(k_omega_Ey_kperp[int((k_omega_Ey_kperp.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto',cmap = plt.get_cmap(cmap))

ax_perp.locator_params(axis='x', nbins=5) #Set number of ticks
ax_perp.tick_params(axis='x',labelsize=labelsize) #Set fontsize of tick labels
ax_perp.xaxis.set_major_locator(MultipleLocator(0.5))

ax_perp.set_xlim(kmin_plot, kmax_plot)
ax_perp.set_ylim(omegamin, omegamax)
ax_perp.set_xlabel("$k_\\bot * d_i$",fontsize=fontsize-1)

ax_perp.plot(a, cfl, linewidth=linewidth, color='black',label='CFL condition')
ax_perp.plot(a,awaves_perp,linewidth=linewidth,color='blue',label='Alfven Speed')
ax_perp.plot(a,fmw_kperp,linewidth=linewidth,color='red',linestyle='--',label='Fast MW')
ax_perp.plot(a,smw_kperp,linewidth=linewidth,color='red',linestyle=':',label='Slow MW')
ax_perp.plot(a,bv_perp,linewidth=linewidth,color='black',linestyle='--',label='V')

ax_perp.text(kmax_plot - 0.3,2.7,'(f)',fontsize=12)

fig.suptitle('$\\Delta r = '+str(int(dx/1000))+'$ km',fontsize=16)

#plt.tight_layout(rect=[-0.01, 0.0, 1, 0.92])
fig.subplots_adjust(wspace=0.05)

#Colorbar
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.855, 0.11, 0.03, 0.77])
fig.colorbar(image, cax=cbar_ax)
cbar_ax.set_ylabel('$E_y^2 (V/m)^2$',fontsize=fontsize-2)


title = 'FFT_'+run
#plt.savefig(path_fig+title+'.eps',dpi=400)
plt.savefig(outputdir+title+'.pdf',dpi=800)
print(outputdir+title+'.pdf')
