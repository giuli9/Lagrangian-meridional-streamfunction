#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import compute_psi_yz as comp
import vars_yz 

if __name__ == "__main__":

    exp="/Users/sauron/Desktop/esperimento_wmed_10days_section1/"
    
    var=vars_yz.Vars(exp)
    jmt_reg=var.jmt_reg
    kmt_reg=var.kmt_reg
    tmask_yz=var.tmask_yz
    
    
    # Costruisco l'array ipb, per tenere traccia dei punti bloccati 
    #   0 1   1 0
    #   1 0   0 1    

    ipb=np.zeros([kmt_reg,jmt_reg])
    for k in range(0,kmt_reg-1):
      for j in range(0,jmt_reg-1):
        if (tmask_yz[k,j] == tmask_yz[k+1,j+1]) and (tmask_yz[k,j+1] == tmask_yz[k+1,j+1]):
          if tmask_yz[k,j]+ tmask_yz[k+1,j] == 1.:
            print("I punti bloccati diagonalmente sono ",k, j)
            ipb[k,j]=1.
                    
    tmask2,j_inisec=comp.section(exp,var)
    mp,mpold=comp.active_points(var,tmask2,ipb)
    psi,ipsi,iref=comp.compute_psi_yz(var,mp,mpold,tmask2,ipb)
#%%					
    yp_reg=var.yp_reg
    zp_reg=var.zp_reg
    yt_reg=var.yt_reg
    zt_reg=var.zt_reg
    y,z=np.meshgrid(yp_reg,-zp_reg)
    fig, ax = plt.subplots(figsize=(15,5))
    vmin = np.nanmin(psi)
    vmax = np.nanmax(psi) 
    levels=np.linspace(-1.05,0.15,17)
    plot=ax.contourf(y,z,psi,cmap="turbo",vmin=vmin,vmax=vmax)
    cbar = plt.colorbar(plot, ax=ax)
    cbar.set_label("Transport (Sv)")
    ticks=[0,-500,-1000,-1500,-2000,-2500,-3000,-4000,-5000]
    ax.set_yticks(ticks)
    z_mask=z.copy()
    z_mask[np.isnan(psi)]=np.nan
    zmin=np.nanmin(z_mask)
    ax.set_ylim([zmin,0])
    #fig.savefig('./Psi.jpg',dpi=500,bbox_inches='tight')
    #cbar.set_ticks([vmin, vmax])  
    #cbar.ax.set_yticklabels([f"{vmin:.3f}", f"{vmax:.3f}"])
    ax.plot([yt_reg[1],yt_reg[1]],[-zt_reg[0],-zt_reg[71]])

	
