#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:41:05 2025

@author: sauron
"""
import numpy as np
import xarray as xr

class Vars:
    def __init__(self, exp):
        # Caricamento dei file NetCDF
        ds = xr.open_dataset(exp + "ariane_statistics_quantitative.nc")
        bs = xr.open_dataset(exp + "mesh_mask_nemo.nc")

        # Lettura degli indici regionali
        self.imt_reg = ds.attrs["imt_reg"]
        self.jmt_reg = ds.attrs["jmt_reg"]
        self.kmt_reg = ds.attrs["kmt_reg"]
        self.imt_reg_start = ds.attrs["imt_reg_start"] - 1
        self.jmt_reg_start = ds.attrs["jmt_reg_start"] - 1
        self.kmt_reg_start = ds.attrs["kmt_reg_start"] - 1
        self.imt_reg_end = ds.attrs["imt_reg_end"]
        self.jmt_reg_end = ds.attrs["jmt_reg_end"]
        self.kmt_reg_end = ds.attrs["kmt_reg_end"]

        # Lettura delle coordinate
        zp = np.squeeze(bs.gdepw_0[:, :])
        yp = np.squeeze(bs.gphif[:, :, 0])

        self.zp_reg = zp[self.kmt_reg_start:self.kmt_reg_end]
        self.yp_reg = yp[self.jmt_reg_start:self.jmt_reg_end]
        
        # lat and lon of T points on C-grid
        zt=np.squeeze(bs.gdept_0[0,:])
        yt=np.squeeze(bs.gphit[0,:,0])
        self.zt_reg=zt[self.kmt_reg_start:self.kmt_reg_end]
        self.yt_reg=yt[self.jmt_reg_start:self.jmt_reg_end]
        
        self.yz_mer = np.squeeze(ds.variables["yz_mer"])
        self.yz_vert = np.squeeze(ds.variables["yz_vert"])
        tmask = ds.tmask

        # Inizializzazione delle maschere
        self.tmask_yz = np.zeros([self.kmt_reg, self.jmt_reg])
        self.pmask_yz = np.ones([self.kmt_reg, self.jmt_reg])  # Valori iniziali a 1
        self.wmask_yz = np.zeros([self.kmt_reg, self.jmt_reg])
        self.vmask_yz = np.zeros([self.kmt_reg, self.jmt_reg])

        # Calcolo della tmask_yz
        for k in range(self.kmt_reg):
            for j in range(self.jmt_reg):
                if np.any(tmask[k, j, :]): 
                    self.tmask_yz[k, j] = 1

        # Calcolo della pmask_yz
        for k in range(self.kmt_reg):
            for j in range(self.jmt_reg):
                if ((self.yz_vert[k, j] == 0 if j < self.jmt_reg else True) 
                    and (self.yz_vert[k, j+1] == 0 if j+1 < self.jmt_reg else True) 
                    and (self.yz_mer[k, j] == 0 if k < self.kmt_reg else True) 
                    and (self.yz_mer[k+1, j] == 0 if k+1 < self.kmt_reg else True)): 
                    self.pmask_yz[k, j] = np.nan

        # Calcolo di wmask_yz
        for k in range(self.kmt_reg - 1):
            for j in range(self.jmt_reg):
                if self.tmask_yz[k, j] + self.tmask_yz[k + 1, j] >= 1:
                    self.wmask_yz[k, j] = 1

        # Calcolo di vmask_yz
        for k in range(self.kmt_reg):
            for j in range(self.jmt_reg - 1):
                if self.tmask_yz[k, j] + self.tmask_yz[k, j + 1] >= 1:
                    self.vmask_yz[k, j] = 1

