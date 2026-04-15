# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 10:02:29 2026

@author: Alejandro
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from PyQt5.QtWidgets import QFileDialog
import re

plt.style.use(["default"])
plt.style.use(["science","notebook"])

try:
    file_path,_ = QFileDialog.getOpenFileName(None,"Open .txt obtained from EXAFS Analysis script")

    df = pd.read_csv(file_path,header = 0,sep='\t')
    
    n = df.keys()
    
    if len(n) == 2:
        print("n = 2")
        
        plt.figure()
        x = df[n[0]]
        y = df[n[1]]
        
        plt.plot(x,y)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Intensity (arb.units)")
        plt.show()
        
    elif len(n) == 8:
        print("n = 8")
        
        expx = df[n[0]]
        expy = df[n[1]]
        
        fx = df[n[2]]
        fy = df[n[3]]

        kchix = df[n[4]]
        kchiy = df[n[5]]

        kfx = df[n[6]]
        kfy = df[n[7]]

        fig, ax = plt.subplots(2, 1, figsize=(10, 8)) 
        
        ax[0].plot(expx,expy, 'black', label="Experimental data")
        ax[0].plot(fx,fy, 'r--', label="Fit")

        ax[0].set_ylabel("|$\\chi(k)$|")
        ax[0].legend(loc='best',frameon=True)
        ax[0].set_xlabel("R ($\\AA$)")
        
        
        ax[1].plot(kchix,kchiy,'b',label='Experimental data')
        ax[1].plot(kfx,kfy,'r--',label='Fit')
        ax[1].set_ylabel("$k^3 \\chi(k)$")
        ax[1].set_xlabel("k ($\\AA^{-1}$)")
        ax[1].legend(frameon=True)
        
        plt.tight_layout()
        plt.show()
        
        
    elif len(n) == 10:
        print("n = 10")
        
        mux = df[n[0]]
        muy = df[n[1]]
        
        pex = df[n[2]]
        pey = df[n[3]]
        
        psx = df[n[4]]
        psy = df[n[5]]
        
        e0x = df[n[6]]
        e0y = df[n[7]]
        
        nx = df[n[8]]
        ny = df[n[9]]
        
        
        fig, ax = plt.subplots(2, 1, figsize=(10, 8)) 
        
        ax[0].plot(mux,muy, 'b', label="Experimental data")
        ax[0].plot(pex,pey, 'g--', label="Pre-edge")
        ax[0].plot(psx,psy, 'r--', label="Post-edge")
        ax[0].axvline(e0x[0],color='k',linestyle='--',label = "$E_0$")
        ax[0].set_ylabel("$\\mu(E)$")
        ax[0].legend(loc='best',frameon=True)
        ax[0].set_xlabel("Energy (eV)")
        
        
        ax[1].plot(nx,ny,'black',label='Normalized data')
        ax[1].set_ylabel("Norm $\\mu(E)$")
        ax[1].set_xlabel("Energy (eV)")
        ax[1].legend(frameon=True)
        
        plt.tight_layout()
        plt.show()
        
    else:
        print("n = 6")
        if re.search("Background",str(df.keys().values)):
            print("Background")
            
            scx = df[n[0]]
            scy = df[n[1]]
            
            bsx = df[n[2]]
            bsy = df[n[3]]    
            
            kchix = df[n[4]]
            kchiy = df[n[5]]
                        
            fig, ax = plt.subplots(2, 1, figsize=(10, 8)) 
            
            ax[0].plot(scx,scy, 'b', label="Experimental data")
            ax[0].plot(bsx,bsy, 'r', label="Background (Spline)")
            ax[0].set_ylabel("$\\mu(E)$")
            ax[0].legend(loc='best',frameon=True)
            ax[0].set_xlabel("Energy (keV)")
            
            
            ax[1].plot(kchix, kchiy,'black',label='Experimental data')
            ax[1].set_ylabel("|$\\chi(k)$|")
            ax[1].set_xlabel("k ($\AA^{-1}$)")
            ax[1].legend(frameon=True)
            
            plt.tight_layout()
            plt.show()
            
            
        elif re.search("Ventana",str(df.keys().values)):
            print("Ventana")
            
            kchix = df[n[0]]
            kchiy = df[n[1]]
            
            vx = df[n[2]]
            vy = df[n[3]]    
            
            chiRX = df[n[4]]
            chiRY = df[n[5]]
            
            fig, ax = plt.subplots(2, 1, figsize=(10, 8)) 
            
            ax[0].plot(kchix, kchiy, 'b', label="Experimental data")
            ax[0].plot(vx, vy, 'r-', label="Window")
            ax[0].set_ylabel("$k^2 \\chi(k)$")
            ax[0].legend(loc='best',frameon=True)
            ax[0].set_xlabel("k ($\AA^{-1}$)")
            
            
            ax[1].plot(chiRX, chiRY, 'black')
            ax[1].set_ylabel("|$\\chi(R)$|")
            ax[1].set_xlabel("R ($\AA$)")
            # ax[1].legend(frameon=True)
            
            plt.tight_layout()
            plt.show()
            
            
        else:
            print("Diferencia")
            
            xmedia_1 = df[n[0]]
            ymedia_1 = df[n[1]]
            
            xmedia_2 = df[n[2]]
            ymedia_2 = df[n[3]]    
            
            diffx = df[n[4]]
            diffy = df[n[5]]
            
            fig, ax = plt.subplots(2, 1, figsize=(8, 6)) 
            
            ax[0].plot(xmedia_1, ymedia_1, 'b', label="1")
            ax[0].plot(xmedia_2, ymedia_2, 'r', label="2")
            ax[0].set_ylabel("$\\mu(E)$")
            ax[0].legend(loc='best',frameon=True)
            
            ax[1].plot(diffx, diffy, 'green', label="Difference (1 - 2)")
            ax[1].set_ylabel("$\Delta$A /-")
            ax[1].set_xlabel("Energy (keV)")
            ax[1].legend(frameon=True)
            
            plt.tight_layout()
            plt.show()
                
    
except:
    print("An exception has occurred")
    
