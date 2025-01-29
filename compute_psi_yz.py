#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import vars_yz 

def section(exp, var):
#     Questa funzione definisce tmask2 come tmask_reg ma con 0 dove sono le sezioni
#     Considera solo un dominio semplice: oltre al coperchio, considera solo sezioni nel piano yz

#     the section indices
    with open(exp+"sections.txt", 'r') as fp:
      nb_sec=len(fp.readlines())

    segind=np.zeros(nb_sec)
    i1=np.zeros(nb_sec,int)
    i2=np.zeros(nb_sec,int)
    j1=np.zeros(nb_sec,int)
    j2=np.zeros(nb_sec,int)
    k1=np.zeros(nb_sec,int)
    k2=np.zeros(nb_sec,int)
    i1_reg=np.zeros(nb_sec,int)
    i2_reg=np.zeros(nb_sec,int)
    j1_reg=np.zeros(nb_sec,int)
    j2_reg=np.zeros(nb_sec,int)
    k1_reg=np.zeros(nb_sec,int)
    k2_reg=np.zeros(nb_sec,int)

    ip=0
    with open(exp+"sections.txt", 'r') as fp:
        for line in fp:
            data = line.split()
            segind[ip]=int(data[0])
            i1[ip]=abs(int(data[1]))-1
            i2[ip]=abs(int(data[2]))-1
            j1[ip]=abs(int(data[3]))-1
            j2[ip]=abs(int(data[4]))-1
            k1[ip]=abs(int(data[5]))-1
            k2[ip]=abs(int(data[6]))-1

            i1_reg[ip]=abs(int(i1[ip]-var.imt_reg_start))
            i2_reg[ip]=abs(int(i2[ip]-var.imt_reg_start))

            j1_reg[ip]=abs(int(j1[ip]-var.jmt_reg_start))
            j2_reg[ip]=abs(int(j2[ip]-var.jmt_reg_start))

            k1_reg[ip]=int(k1[ip]-var.kmt_reg_start)
            k2_reg[ip]=int(k2[ip]-var.kmt_reg_start)
            ip+=1
                    
    tmask2 = var.tmask_yz.copy()        
    for ii in range(0,nb_sec):
        if ii==0:
            j_inisec=np.array([j1_reg[ii],j2_reg[ii]])
        if segind[ii] > 0:
            if j1_reg[ii] == j2_reg[ii]:
                for k in range(k1_reg[ii],k2_reg[ii]+1):
                    if var.tmask_yz[k,j1_reg[ii]] == 1.:
                        tmask2[k,j1_reg[ii]] = 0.
                    
            #elif (segind[ii] > 0) & (j1_reg[ii] == j2_reg[ii]):
            #    for i in range(i1_reg[ii],i2_reg[ii]+1):
            #        if tmask_yz[k,j1_reg[ii]] == 1:
            #            tmask2[k,j1_reg[ii]] = 0
            elif (k1_reg[ii] == -1) and (k2_reg[ii] == k1_reg[ii]):
                print("Section {} is a lid".format(int(segind[ii])))
                continue    
  
    return tmask2,j_inisec
                    
def active_points(var, tmask2, ipb):
# Questa funzione trova i punti attivi kp0 e jp0, e costruisce in base a questo
# l'array mp, che ha 1 in tutti i punti attivi e 0 altrimenti.
# Per calcolare mp, ho bisogno di ipb e tmask2      
# return mp, mpold
# mpold è il numero di tutti i punti attivi
    jmt_reg=var.jmt_reg
    kmt_reg=var.kmt_reg
    yz_mer=var.yz_mer
    yz_vert=var.yz_vert
    
    div=np.zeros([kmt_reg,jmt_reg])
    for k in range(1,var.kmt_reg):
        for j in range(1,jmt_reg):
            div[k,j]=yz_vert[k,j]-yz_vert[k-1,j]+yz_mer[k,j]-yz_mer[k,j-1]
            
    kp0=0
    jp0=0
    divmax=np.max(abs(div))        
    for k in range(1,kmt_reg-1):
        for j in range(1,jmt_reg-1):
            if div[k,j]<=divmax and div[k+1,j]<=divmax and div[k,j+1]<=divmax and div[k+1,j+1]<=divmax:
                if yz_mer[k,j]*yz_mer[k+1,j]!=0 or yz_vert[k,j]*yz_vert[k,j+1]!=0:
                    kp0=k
                    jp0=j
                    break
                
        if kp0!=0 and jp0!=0:
            print("Active points found: z_indx=",kp0," y_indx=",jp0) 
            break
        
    mp=np.zeros([kmt_reg,jmt_reg]) 
    mp[kp0,jp0]=1.

    print('\n')
    print('Computing psi mask (mp)');
    mpold=0
    totmp=1
    while totmp > mpold: 
        print('\n')
        print('Totmp: ',totmp,' mpold: ',mpold);
        mpold=totmp
        
        for k in range(0,kmt_reg-1):
            for j in range(0,jmt_reg-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k+1,j]==0. and (tmask2[k+1,j]==1. or tmask2[k+1,j+1]==1. ):
                        mp[k+1,j]=1. 
                    if mp[k,j+1]==0. and (tmask2[k,j+1]==1. or tmask2[k+1,j+1]==1.):
                        mp[k,j+1]=1.
                        
        for k in range(1,kmt_reg-1):
            for j in range(1,jmt_reg-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k-1,j]==0. and (tmask2[k,j]==1. or tmask2[k,j+1]==1. ):
                        mp[k-1,j]=1. 
                    if mp[k,j-1]==0. and (tmask2[k,j]==1. or tmask2[k+1,j]==1.):
                        mp[k,j-1]=1.
                        
        for k in range(0,kmt_reg-1):
            for j in range(jmt_reg-2,-1,-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k+1,j]==0. and (tmask2[k+1,j]==1. or tmask2[k+1,j+1]==1. ):
                        mp[k+1,j]=1. 
                    if mp[k,j+1]==0. and (tmask2[k,j+1]==1. or tmask2[k+1,j+1]==1.):
                        mp[k,j+1]=1.
                        
        for k in range(1,kmt_reg-1):
            for j in range(jmt_reg-2,0,-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k-1,j]==0. and (tmask2[k,j]==1. or tmask2[k,j+1]==1. ):
                        mp[k-1,j]=1. 
                    if mp[k,j-1]==0. and (tmask2[k,j]==1. or tmask2[k+1,j]==1.):
                        mp[k,j-1]=1.
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(0,jmt_reg-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k+1,j]==0. and (tmask2[k+1,j]==1. or tmask2[k+1,j+1]==1. ):
                        mp[k+1,j]=1. 
                    if mp[k,j+1]==0. and (tmask2[k,j+1]==1. or tmask2[k+1,j+1]==1.):
                        mp[k,j+1]=1.
                        
        for k in range(kmt_reg-2,0,-1):
            for j in range(1,jmt_reg-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k-1,j]==0. and (tmask2[k,j]==1. or tmask2[k,j+1]==1. ):
                        mp[k-1,j]=1. 
                    if mp[k,j-1]==0. and (tmask2[k,j]==1. or tmask2[k+1,j]==1.):
                        mp[k,j-1]=1.
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(jmt_reg-2,-1,-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k+1,j]==0. and (tmask2[k+1,j]==1. or tmask2[k+1,j+1]==1. ):
                        mp[k+1,j]=1. 
                    if mp[k,j+1]==0. and (tmask2[k,j+1]==1. or tmask2[k+1,j+1]==1.):
                        mp[k,j+1]=1.
                        
        for k in range(kmt_reg-2,0,-1):
            for j in range(jmt_reg-2,0,-1):
                if ipb[k,j]==0. and mp[k,j]==1.:
                    if mp[k-1,j]==0. and (tmask2[k,j]==1. or tmask2[k,j+1]==1. ):
                        mp[k-1,j]=1. 
                    if mp[k,j-1]==0. and (tmask2[k,j]==1. or tmask2[k+1,j]==1.):
                        mp[k,j-1]=1.
                        
        
        totmp=sum(sum(mp)) #Continuo finchè non trovo tutti i punti attivi
    return mp, mpold

def compute_psi_yz(var,mp, mpold, tmask2, ipb):
# Prima trovo il punto da cui parte l'integrazione e in cui psi[kref_psi,jref_psi]=0
# poi costruisco iref che mi dice quanti e quali punto sono connessi a quello di referenza,
#    e quindi su quali punto posso integrare
    jmt_reg=var.jmt_reg
    kmt_reg=var.kmt_reg
    yz_mer=var.yz_mer
    yz_vert=var.yz_vert
    pmask_yz=var.pmask_yz
    wmask_yz=var.wmask_yz
    vmask_yz=var.vmask_yz
    
    psi=np.zeros([kmt_reg,jmt_reg])-1.e12
    ipsi=np.zeros([kmt_reg,jmt_reg])
    
    iref=np.zeros([kmt_reg,jmt_reg])
    
    #Selezioni gli indici che definiscono il punto dal quale parto a integrare la psi
    kref_psi=0
    jref_psi=1
    iref[kref_psi,jref_psi]=1
    
    irefold=0
    totiref=1
    
    while totiref != irefold: 
        
        print("Punti attivi: ",irefold," per un totale di ",totiref)
        irefold=totiref
        
        for k in range(0,kmt_reg-1):
            for j in range(0,jmt_reg-1):
                if iref[k,j]==1.:
                    if iref[k+1,j]==0. and tmask2[k+1,j]==1. and tmask2[k+1,j+1]==1.:
                        iref[k+1,j]=1. 
                    if iref[k,j+1]==0. and tmask2[k,j+1]==1. and tmask2[k+1,j+1]==1.:
                        iref[k,j+1]=1.
                        
        for k in range(1,kmt_reg-1):
            for j in range(1,jmt_reg-1):
                if iref[k,j]==1.:
                    if iref[k-1,j]==0. and tmask2[k,j]==1. and tmask2[k,j+1]==1.:
                        iref[k-1,j]=1. 
                    if iref[k,j-1]==0. and tmask2[k,j]==1. and tmask2[k+1,j]==1.:
                        iref[k,j-1]=1.
                        
        for k in range(0,kmt_reg-1):
            for j in range(jmt_reg-2,-1,-1):
                if iref[k,j]==1.:
                    if iref[k+1,j]==0. and tmask2[k+1,j]==1. and tmask2[k+1,j+1]==1.:
                        iref[k+1,j]=1. 
                    if iref[k,j+1]==0. and tmask2[k,j+1]==1. and tmask2[k+1,j+1]==1.:
                        iref[k,j+1]=1.
                        
        for k in range(1,kmt_reg-1):
            for j in range(jmt_reg-2,0,-1):
                if iref[k,j]==1.:
                    if iref[k-1,j]==0. and tmask2[k,j]==1. and tmask2[k,j+1]==1.:
                        iref[k-1,j]=1. 
                    if iref[k,j-1]==0. and tmask2[k,j]==1. and tmask2[k+1,j]==1.:
                        iref[k,j-1]=1.
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(0,jmt_reg-1):
                if iref[k,j]==1.:
                    if iref[k+1,j]==0. and tmask2[k+1,j]==1. and tmask2[k+1,j+1]==1.:
                        iref[k+1,j]=1. 
                    if iref[k,j+1]==0. and tmask2[k,j+1]==1. and tmask2[k+1,j+1]==1.:
                        iref[k,j+1]=1.
                        
        for k in range(kmt_reg-2,0,-1):
            for j in range(1,jmt_reg-1):
                if iref[k,j]==1.:
                    if iref[k-1,j]==0. and tmask2[k,j]==1. and tmask2[k,j+1]==1.:
                        iref[k-1,j]=1. 
                    if iref[k,j-1]==0. and tmask2[k,j]==1. and tmask2[k+1,j]==1.:
                        iref[k,j-1]=1.
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(jmt_reg-2,-1,-1):
                if iref[k,j]==1.:
                    if iref[k+1,j]==0. and tmask2[k+1,j]==1. and tmask2[k+1,j+1]==1.:
                        iref[k+1,j]=1. 
                    if iref[k,j+1]==0. and tmask2[k,j+1]==1. and tmask2[k+1,j+1]==1.:
                        iref[k,j+1]=1.
                        
        for k in range(kmt_reg-2,0,-1):
            for j in range(jmt_reg-2,0,-1):
                if iref[k,j]==1.:
                    if iref[k-1,j]==0. and tmask2[k,j]==1. and tmask2[k,j+1]==1.:
                        iref[k-1,j]=1. 
                    if iref[k,j-1]==0. and tmask2[k,j]==1. and tmask2[k+1,j]==1.:
                        iref[k,j-1]=1.
                        
        
        totiref=sum(sum(iref)) #Continuo finchè non trovo tutti i punti attivi
    
    kp1=56 #41
    jp1=48 #54
    for k in range(0,kmt_reg):
        for j in range(0,jmt_reg):
            if iref[k,j]==1. and mp[k,j]==1. and kp1==0:
                kp1=k
                jp1=j
                
    print("L'integrazione inizia agli indici k=", kp1," e y=",jp1)
    print("\n")
    print("qui yz_vert=",yz_vert[kp1,jp1].values,"\n")
    print("qui yz_mer=",yz_mer[kp1,jp1].values,"\n")
    
    
# Definisco psi e ipsi. Quest'ultima tiene traccia dei punti
# del dominio già integrati

    psi[kp1,jp1]=0.
    ipsi[kp1,jp1]=1.
    totipsi=1
    
    while totipsi<mpold:
        print("Totipsi è ", totipsi,"\n")
        print("Mpold è ", mpold, "\n")
        
        for k in range(0,kmt_reg):
            for j in range(0,jmt_reg-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j+1]==0 and mp[k,j+1]==1. and wmask_yz[k,j+1]==1.:
                        psi[k,j+1]=psi[k,j]+yz_vert[k,j+1]
                        ipsi[k,j+1]=1. 
                        
        for k in range(0,kmt_reg):
            for j in range(1,jmt_reg):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j-1]==0. and mp[k,j-1]==1. and wmask_yz[k,j]==1.:
                        psi[k,j-1]=psi[k,j]-yz_vert[k,j]
                        ipsi[k,j-1]=1.
        
        for k in range(0,kmt_reg-1):
            for j in range(jmt_reg-1,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k+1,j]==0. and mp[k+1,j]==1. and vmask_yz[k+1,j]==1.:
                        psi[k+1,j]=psi[k,j]-yz_mer[k+1,j]
                        ipsi[k+1,j]=1.
                        
        for k in range(1,kmt_reg):
            for j in range(jmt_reg-1,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k-1,j]==0. and mp[k-1,j]==1. and vmask_yz[k,j]==1.:
                        psi[k-1,j]=psi[k,j]+yz_mer[k,j]
                        ipsi[k-1,j]=1.
                        
        for k in range(0,kmt_reg):
            for j in range(jmt_reg-2,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j+1]==0 and mp[k,j+1]==1. and wmask_yz[k,j+1]==1.:
                        psi[k,j+1]=psi[k,j]+yz_vert[k,j+1]
                        ipsi[k,j+1]=1. 
                        
        for k in range(0,kmt_reg):
            for j in range(jmt_reg-1,0,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j-1]==0. and mp[k,j-1]==1. and wmask_yz[k,j]==1.:
                        psi[k,j-1]=psi[k,j]-yz_vert[k,j]
                        ipsi[k,j-1]=1.        
                        
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(jmt_reg-1,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k+1,j]==0. and mp[k+1,j]==1. and vmask_yz[k+1,j]==1.: 
                        psi[k+1,j]=psi[k,j]-yz_mer[k+1,j] 
                        ipsi[k+1,j]=1.
                        
        for k in range(kmt_reg-1,0,-1):
            for j in range(jmt_reg-1,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k-1,j]==0. and mp[k-1,j]==1. and vmask_yz[k,j]==1.:
                        psi[k-1,j]=psi[k,j]+yz_mer[k,j]
                        ipsi[k-1,j]=1.    
                        
        for k in range(kmt_reg-1,-1,-1):
            for j in range(jmt_reg-2,-1,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j+1]==0 and mp[k,j+1]==1. and wmask_yz[k,j+1]==1.:
                        psi[k,j+1]=psi[k,j]+yz_vert[k,j+1]
                        ipsi[k,j+1]=1. 
                        
        for k in range(kmt_reg-1,-1,-1):
            for j in range(jmt_reg-1,0,-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j-1]==0. and mp[k,j-1]==1. and wmask_yz[k,j]==1.:
                        psi[k,j-1]=psi[k,j]-yz_vert[k,j]
                        ipsi[k,j-1]=1.
                        
                        
        for k in range(kmt_reg-2,-1,-1):
            for j in range(0,jmt_reg):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k+1,j]==0. and mp[k+1,j]==1. and vmask_yz[k+1,j]==1.: 
                        psi[k+1,j]=psi[k,j]-yz_mer[k+1,j] 
                        ipsi[k+1,j]=1.
                        
        for k in range(kmt_reg-1,0,-1):
            for j in range(0,jmt_reg):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k-1,j]==0. and mp[k-1,j]==1. and vmask_yz[k,j]==1.:
                        psi[k-1,j]=psi[k,j]+yz_mer[k,j]
                        ipsi[k-1,j]=1.    
                        
        for k in range(kmt_reg-1,-1,-1):
            for j in range(0,jmt_reg-1):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j+1]==0 and mp[k,j+1]==1. and wmask_yz[k,j+1]==1.:
                        psi[k,j+1]=psi[k,j]+yz_vert[k,j+1]
                        ipsi[k,j+1]=1. 
                        
        for k in range(kmt_reg-1,-1,-1):
            for j in range(1,jmt_reg):
                if ipsi[k,j]==1. and ipb[k,j]==0.:
                    if ipsi[k,j-1]==0. and mp[k,j-1]==1. and wmask_yz[k,j]==1.:
                        psi[k,j-1]=psi[k,j]-yz_vert[k,j]
                        ipsi[k,j-1]=1.
                        
                        
        totipsi=sum(sum(ipsi))
        
    psi=psi/1.e6; psi[psi==-1.e6]=np.nan; psi=psi*pmask_yz
    
    
    return psi,ipsi, iref
                        
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
                    
    tmask2,j_inisec=section(exp,var)
    mp,mpold=active_points(var,tmask2)
    psi,ipsi,iref=compute_psi_yz(var,mp,mpold,tmask2)

    