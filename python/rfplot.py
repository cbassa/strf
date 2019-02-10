#!/usr/bin/env python3
import numpy as np
import ppgplot
from strf.rfio import Spectrogram
from strf.cmap import cm_cool, cm_viridis, cm_heat

if __name__ == "__main__":
    # Settings
    path = "/data2/satobs/radio/graves_dt/bin"
    prefix = "2018-12-31T09:34:11"
    ifile = 48000
    nsub = 900
    path = "/data2/satobs/radio/20190207"
    prefix = "2019-02-07T08:26:59"
    ifile = 0
    nsub = 21600

    # Read spectrogram
    s = Spectrogram.read_bin_files(path, prefix, ifile, nsub, 4171)
    
    zmin = 0.0
    zmax = 0.00008
    zmax = 0.02


    # Start ppgplot
    ppgplot.pgopen("/xs")
    ppgplot.pgsch(0.8)
    ppgplot.pgask(0)

    # Default settings
    xmin, ymin = 0, 0
    xmax, ymax = s.nsub, s.nchan
    mode, posn, click = 0, 0, 0
    x0, y0 = 0.0, 0.0
    cmap = 2
    
    # Image limits
    tr = np.array([-0.5, 1.0, 0.0, -0.5, 0.0, 1.0])
    
    # Forever loop
    redraw = True
    while True:
        # Redraw
        if redraw:
            # Set up image view port
            ppgplot.pgpage()
            ppgplot.pgsci(1)
            ppgplot.pgsvp(0.1, 0.95, 0.1, 0.95)
            ppgplot.pgswin(xmin ,xmax, ymin, ymax)

            # Plot spectrogram
            if cmap==3:
                ppgplot.pggray(s.z, s.nsub, s.nchan, 0, s.nsub-1, 0, s.nchan-1, zmin, zmax, tr)
            else:
                if cmap==0:
                    ppgplot.pgctab(cm_cool[0], cm_cool[1], cm_cool[2], cm_cool[3], len(cm_cool[0]), 1.0, 0.5)
                elif cmap==1: 
                    ppgplot.pgctab(cm_heat[0], cm_heat[1], cm_heat[2], cm_heat[3], len(cm_heat[0]), 1.0, 0.5)
                elif cmap==2:
                    ppgplot.pgctab(cm_viridis[0], cm_viridis[1], cm_viridis[2], cm_viridis[3], len(cm_viridis[0]), 1.0, 0.5)
                ppgplot.pgimag(s.z, s.nsub, s.nchan, 0, s.nsub-1, 0, s.nchan-1, zmax, zmin, tr)

            # Plot pixel axis
            ppgplot.pgbox("CTSM1", 0., 0, "CTSM1", 0., 0)

            # Time axis
            ppgplot.pgbox("B", 0., 0, "", 0., 0);

            # Frequency axis
            fmin, fmax = np.min(s.freq)*1e-6, np.max(s.freq)*1e-6
            fcen = np.round(500.0*(fmax+fmin))/1000.0 # Round to nearest kHz
            fmin -= fcen
            fmax -= fcen
            ppgplot.pgswin(xmin, xmax, fmin, fmax)
            ppgplot.pgbox("", 0., 0, "BCTSN", 0., 0)
            ppgplot.pglab("UT Date", "Frequency - {:.3f} MHz".format(fcen), "")

            # Back to pixels
            ppgplot.pgswin(xmin, xmax, ymin, ymax)
            
            
            # Reset flag
            redraw = False

        # Get cursor
        x, y, char = ppgplot.pgband(mode, posn, x0, y0)
        c = char.decode("utf-8")
        
        # Quit
        if c=="q":
            break

        # Toggle colormap
        if c=="C":
            cmap+=1
            if cmap>3:
                cmap=1
            redraw = True
            continue
        
        # Reset
        if c=="r":
            xmin, ymin = 0.0, 0.0
            xmax, ymax = s.nsub, s.nchan
            redraw = True
            continue

        # Zoom
        if c=="z":
            click, mode = 1, 2

        # Execute zoom
        if c=="A":
            if click==0:
                click = 1
            elif (click==1) and (mode==2):
                xmin, xmax = min(x0, x), max(x0, x)
                ymin, ymax = min(y0, y), max(y0, y)
                click, mode = 0, 0
                redraw = True
            else:
                click, mode = 0, 0
                redraw = True
                
        # Save past cursor position
        x0, y0 = x, y
                
        
            
    ppgplot.pgend()
