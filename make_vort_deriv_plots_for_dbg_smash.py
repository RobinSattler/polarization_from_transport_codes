import fileinput
import math
import numpy as np
import sys
import os
import pickle
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm
hbarc=0.197326


time_min=15.00
time_max=15.00

xmax=15
ymax=15
zmax=15

minlogval=-3 #minimum value to display in log plots

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files!=3):
   print ('Syntax: ./make_cg_plots.py <pickle binary data file> <output directory> <plot common title>')
   sys.exit(1)

inputfile=sys.argv[1]
od=sys.argv[2]
common_title=sys.argv[3]
if(not os.path.exists(od)):
  os.mkdir(od)

with open(inputfile,"rb") as pi:
     data=pickle.load(pi)

tt,xx,yy,zz,vx,vy,vz,temp,bt,bx,by,bz,dbt_dx,dbt_dy,dbt_dz,dbx_dt,dby_dt,dbz_dt,dbx_dy,dbx_dz,dby_dx,dby_dz,dbz_dx,dbz_dy=data[:]

omega_tx=0.5*hbarc*(dbt_dx-dbx_dt)
omega_ty=0.5*hbarc*(dbt_dy-dby_dt)
omega_tz=0.5*hbarc*(dbt_dz-dbz_dt)

omega_yz=0.5*hbarc*(dby_dz-dbz_dy)
omega_zx=0.5*hbarc*(dbz_dx-dbx_dz)
omega_xy=0.5*hbarc*(dbx_dy-dby_dx)

#print(str(xx))
#lx=len(xx)
#print(str(lx))
#ic=int(math.floor(lx/2))
#print(str(ic))
#print("x dir: "+str(vx[34,ic-1:ic+2,ic,ic]))
#print("y dir: "+str(vy[34,ic,ic-1:ic+2,ic]))
#print("z dir: "+str(vz[34,ic,ic,ic-1:ic+2]))
#sys.exit(0)

im=np.argmin(abs(xx+xmax))
jm=np.argmin(abs(yy+ymax))
km=np.argmin(abs(zz+zmax))
i0=np.argmin(abs(xx))
j0=np.argmin(abs(yy))
k0=np.argmin(abs(zz))
ip=np.argmin(abs(xx-xmax))+1
jp=np.argmin(abs(yy-ymax))+1
kp=np.argmin(abs(zz-zmax))+1

#print(str(xx[im])+"  "+str(yy[jm])+"  "+str(zz[km]))
#print(str(xx[i0])+"  "+str(yy[j0])+"  "+str(zz[k0]))
#print(str(xx[ip])+"  "+str(yy[jp])+"  "+str(zz[kp]))
#sys.exit(0)

dx=xx[1]-xx[0]
dx2=dx/2
dy=yy[1]-yy[0]
dy2=dy/2
dz=zz[1]-zz[0]
dz2=dz/2

XYE=[xx[im]-dx2, xx[ip-1]+dx2, yy[jm]-dy2, yy[jp-1]+dy2]
XZE=[xx[im]-dx2, xx[ip-1]+dx2, zz[km]-dz2, zz[kp-1]+dz2]
YZE=[yy[jm]-dy2, yy[jp-1]+dy2, zz[km]-dz2, zz[kp-1]+dz2]

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 18
fig_size[1] = 8
plt.rcParams["figure.figsize"]=fig_size


for it in range(len(tt)):
  if(tt[it]<time_min):
    continue
  if(tt[it]>time_max):
    sys.exit(0)
  print("*****\n\nDoing timestep "+str(it)+", t="+'{:4.2f}'.format(tt[it]))


  plt.suptitle(common_title+", "+"$v_x$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=vx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=vx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(vx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^x(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=vx[it,im:ip,j0,km:kp].max()
  vmin_tmp=vx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(vx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^x(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=vx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=vx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(vx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^x(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/vx_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')



  plt.suptitle(common_title+", "+"$v_y$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=vy[it,im:ip,jm:jp,k0].min()
  vmax_tmp=vy[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(vy[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^y(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=vy[it,im:ip,j0,km:kp].max()
  vmin_tmp=vy[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(vy[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^y(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=vy[it,i0,jm:jp,km:kp].max()
  vmin_tmp=vy[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(vy[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^y(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.tight_layout()


  plt.savefig(od+"/vy_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
    

  plt.suptitle(common_title+", "+"$v_z$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=vz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=vz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(vz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=vz[it,im:ip,j0,km:kp].max()
  vmin_tmp=vz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(vz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=vz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=vz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(vz[it,i0,jm:jp,km:kp].transpose(),extent=YZE, origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$v^z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='c units',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/vz_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
    

  plt.suptitle(common_title+", "+"$T$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  plt.imshow(1000*temp[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$T(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='MeV',cax=cax,format="%6.3f")

  plt.subplot(132)
  plt.imshow(1000*temp[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$T(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='MeV',cax=cax,format="%6.3f")

  plt.subplot(133)
  plt.imshow(1000*temp[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$T(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='MeV',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/temp_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
    
  plt.suptitle(common_title+", "+r"$\beta_t$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=bt[it,im:ip,jm:jp,k0].min()
  vmax_tmp=bt[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bt[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_t(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=bt[it,im:ip,j0,km:kp].max()
  vmin_tmp=bt[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bt[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_t(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=bt[it,i0,jm:jp,km:kp].max()
  vmin_tmp=bt[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bt[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_t(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/bt_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
    
  plt.suptitle(common_title+", "+r"$\beta_x$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=bx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=bx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_x(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=bx[it,im:ip,j0,km:kp].max()
  vmin_tmp=bx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_x(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=bx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=bx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_x(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/bx_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')


  plt.suptitle(common_title+", "+r"$\beta_y$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=by[it,im:ip,jm:jp,k0].min()
  vmax_tmp=by[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(by[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_y(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=by[it,im:ip,j0,km:kp].max()
  vmin_tmp=by[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(by[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_y(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=by[it,i0,jm:jp,km:kp].max()
  vmin_tmp=by[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(by[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_y(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/by_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')


  plt.suptitle(common_title+", "+r"$\beta_z$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(231)
  vmin_tmp=bz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=bz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(bz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(232)
  vmax_tmp=bz[it,im:ip,j0,km:kp].max()
  vmin_tmp=bz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(bz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(233)
  vmax_tmp=bz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=bz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(bz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\beta_z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(234)
  plt.imshow(np.log10(np.abs(bz[it,im:ip,jm:jp,k0].transpose())+1.e-10),vmin=minlogval,extent=XYE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$Log_{10}\beta_z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(235)
  plt.imshow(np.log10(np.abs(bz[it,im:ip,j0,km:kp].transpose())+1.e-10),vmin=minlogval,extent=XZE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$Log_{10}\beta_z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(236)
  plt.imshow(np.log10(np.abs(bz[it,i0,jm:jp,km:kp].transpose())+1.e-10),vmin=minlogval,extent=YZE,origin='lower',cmap='gnuplot',aspect='equal',interpolation=None)
  plt.title(r"$Log_{10}\beta_z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout(pad=3)

  plt.savefig(od+"/bz_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')



  plt.suptitle(common_title+", "+r"$\partial \beta_t/\partial x$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbt_dx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbt_dx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial x(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbt_dx[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbt_dx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial x(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbt_dx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbt_dx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial x(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbt_x_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_t/\partial y$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbt_dy[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbt_dy[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dy[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial y(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbt_dy[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbt_dy[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dy[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial y(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbt_dy[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbt_dy[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dy[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial y(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbt_y_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_t/\partial z$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbt_dz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbt_dz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbt_dz[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbt_dz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbt_dz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbt_dz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbt_dz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_t/\partial z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbt_z_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_x/\partial t$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbx_dt[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbx_dt[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dt[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial t(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbx_dt[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbx_dt[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dt[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial t(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbx_dt[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbx_dt[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dt[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial t(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbx_t_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_y/\partial t$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dby_dt[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dby_dt[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dt[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial t(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dby_dt[it,im:ip,j0,km:kp].max()
  vmin_tmp=dby_dt[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dt[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial t(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dby_dt[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dby_dt[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dt[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial t(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dby_t_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_z/\partial t$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbz_dt[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbz_dt[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbz_dt[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial t(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbz_dt[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbz_dt[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbz_dt[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial t(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbz_dt[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbz_dt[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbz_dt[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial t(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbz_t_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_x/\partial y$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbx_dy[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbx_dy[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dy[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial y(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbx_dy[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbx_dy[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dy[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial y(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbx_dy[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbx_dy[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dy[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial y(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbx_y_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_x/\partial z$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbx_dz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbx_dz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbx_dz[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbx_dz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbx_dz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbx_dz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbx_dz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_x/\partial z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbx_z_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_y/\partial z$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dby_dz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dby_dz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial z(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dby_dz[it,im:ip,j0,km:kp].max()
  vmin_tmp=dby_dz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial z(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dby_dz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dby_dz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial z(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dby_z_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_y/\partial x$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dby_dx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dby_dx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial x(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dby_dx[it,im:ip,j0,km:kp].max()
  vmin_tmp=dby_dx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial x(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dby_dx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dby_dx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dby_dx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_y/\partial x(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dby_x_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_z/\partial x$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbz_dx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbz_dx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(dbz_dx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial x(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbz_dx[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbz_dx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(dbz_dx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial x(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbz_dx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbz_dx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(dbz_dx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial x(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbz_x_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\partial \beta_z/\partial y$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=dbz_dy[it,im:ip,jm:jp,k0].min()
  vmax_tmp=dbz_dy[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(dbz_dy[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial y(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=dbz_dy[it,im:ip,j0,km:kp].max()
  vmin_tmp=dbz_dy[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(dbz_dy[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial y(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=dbz_dy[it,i0,jm:jp,km:kp].max()
  vmin_tmp=dbz_dy[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(dbz_dy[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\partial \beta_z/\partial y(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/dbz_y_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\omega_{tx}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_tx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_tx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_tx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tx}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_tx[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_tx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_tx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tx}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_tx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_tx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_tx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tx}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_tx_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')


  plt.suptitle(common_title+", "+r"$\omega_{ty}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_ty[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_ty[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_ty[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{ty}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_ty[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_ty[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_ty[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{ty}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_ty[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_ty[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_ty[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{ty}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_ty_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\omega_{tz}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_tz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_tz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_tz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tz}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_tz[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_tz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_tz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tz}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_tz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_tz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_tz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{tz}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_tz_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')

  plt.suptitle(common_title+", "+r"$\omega_{yz}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_yz[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_yz[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_yz[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{yz}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_yz[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_yz[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_yz[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{yz}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_yz[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_yz[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_yz[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{yz}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_yz_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')


  plt.suptitle(common_title+", "+r"$\omega_{zx}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_zx[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_zx[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_zx[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{zx}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_zx[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_zx[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_zx[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{zx}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_zx[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_zx[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_zx[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{zx}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_zx_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
  

  plt.suptitle(common_title+", "+r"$\omega_{xy}$"+", t="+'{:4.2f}'.format(tt[it]))

  plt.subplot(131)
  vmin_tmp=omega_xy[it,im:ip,jm:jp,k0].min()
  vmax_tmp=omega_xy[it,im:ip,jm:jp,k0].max()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)

  plt.imshow(omega_xy[it,im:ip,jm:jp,k0].transpose(),extent=XYE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{xy}(x,y,0)$")
  plt.xlabel('x [fm]')
  plt.ylabel('y [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(132)
  vmax_tmp=omega_xy[it,im:ip,j0,km:kp].max()
  vmin_tmp=omega_xy[it,im:ip,j0,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2. 
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_xy[it,im:ip,j0,km:kp].transpose(),extent=XZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{xy}(x,0,z)$")
  plt.xlabel('x [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.subplot(133)
  vmax_tmp=omega_xy[it,i0,jm:jp,km:kp].max()
  vmin_tmp=omega_xy[it,i0,jm:jp,km:kp].min()
  if((vmax_tmp > 0) and (vmin_tmp <0)):
      vcenter_tmp=0.
      cmap_tmp='seismic'
      if(vmax_tmp > -vmin_tmp):
          vmin_tmp=-vmax_tmp
      else:
          vmax_tmp=-vmin_tmp
  else:
      cmap_tmp='gnuplot'
      vcenter_tmp=(vmax_tmp+vmin_tmp)/2.
  norm = DivergingNorm(vmin=vmin_tmp, vcenter=vcenter_tmp, vmax=vmax_tmp)
  plt.imshow(omega_xy[it,i0,jm:jp,km:kp].transpose(),extent=YZE,origin='lower',cmap=cmap_tmp,aspect='equal',norm=norm,interpolation=None)
  plt.title(r"$\omega_{xy}(0,y,z)$")
  plt.xlabel('y [fm]')
  plt.ylabel('z [fm]')
  ax=plt.gca()
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(label='$fm^{-1}GeV^{-1}$',cax=cax,format="%6.3f")

  plt.tight_layout()

  plt.savefig(od+"/omega_xy_time_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
  plt.close('all')
