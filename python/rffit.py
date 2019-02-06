#!/usr/bin/env python
import ephem
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime
from scipy import interpolate

class satellite:
    """Satellite class"""

    def __init__(self, tle0, tle1, tle2):
        """Define a satellite"""

        self.tle0 = tle0
        self.tle1 = tle1
        self.tle2 = tle2
        if tle0[:2]=="0 ":
            self.name = tle0[2:]
        else:
            self.name = tle0
        self.id = int(tle1.split(" ")[1][:5])

def compute_altitude_and_velocity(observer, satellite, times):
    alt = []
    v = []
    for t in times:
        observer.date = t.datetime
        satellite.compute(observer)
        alt.append(satellite.alt)
        v.append(satellite.range_velocity)

    return np.array(alt), np.array(v)

if __name__ == "__main__":
    # Settings
    c = 299792458.0

    # Set RX
    obsrx = ephem.Observer()
    obsrx.lat = "52.8119"
    obsrx.lon = "6.3961"
    obsrx.elevation = 10

    # Set TX
    obstx = ephem.Observer()
    obstx.lat = "47.348"
    obstx.lon = "5.5151"
    obstx.elevation = 100

    # Read data
    d = np.loadtxt("data/2019-02-03T22:18:02_0143.050.dat")
    mjd = d[:, 0]
    freq = d[:, 1]
    times = Time(mjd, format="mjd", scale="utc")

    tmin, tmax = np.min(times), np.max(times)
    timesmod = Time(np.linspace(tmin.mjd, tmax.mjd, 90), format="mjd", scale="utc")
    mjdmod = timesmod.mjd
    
    # Read satellites
    with open("bulk.tle", "r") as f:
        lines = f.readlines()
        satellites = [satellite(lines[i], lines[i+1], lines[i+2])
                      for i in range(0, len(lines), 3)]

    # Read frequencies
    freqs = dict()
    with open("frequencies.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            norad_cat_id = int(line.split("  ")[0].encode("ascii"))
            freqs[norad_cat_id] = float(line.split("  ")[1])*1e6

    # Axis means
    freq0 = 143050000
    mjd0 = np.floor(np.min(mjd))

    # Generate x and y positions
    x = 86400.0*(mjd-mjd0)
    y = freq-freq0
            
    plt.figure(figsize=(20, 10))
    
    # Loop over satellites
    for sat in satellites:
        # Skip satellites not in the list of frequencies
        if sat.id not in freqs:
            continue

        # Set satellite
        satellite = ephem.readtle(sat.tle0, sat.tle1, sat.tle2)

        # See if object is above the horizon
        alt, v = compute_altitude_and_velocity(obstx, satellite, [tmin, tmax])
        
        # Compute velocities
        if alt[0]>0.0 or alt[1]>0.0:
            altrx, vrx = compute_altitude_and_velocity(obsrx, satellite, timesmod)
            alttx, vtx = compute_altitude_and_velocity(obstx, satellite, timesmod)
            freqmod = freqs[sat.id]*(1.0-vrx/c)*(1.0-vtx/c)

            # Cmopute offsets
            xm = 86400.0*(mjdmod-mjd0)
            ym = freqmod-freq0
            zm = alttx
            
            # Define interpolating function
            ffreq = interpolate.interp1d(xm, ym)
            falt = interpolate.interp1d(xm, zm)
            dy = y-ffreq(x)
            cmatch = (np.abs(dy)<30) & (falt(x)>0.0) & (np.abs(y)>100.0)
            if np.sum(cmatch)>0:
                print("%05d %6d %8.3f %8.3f"%(sat.id, np.sum(cmatch), np.mean(dy[cmatch]), np.std(dy[cmatch])))
                plt.plot(x[cmatch], y[cmatch], "ro")
            plt.plot(xm, ym, alpha=0.1, color="C0")

    plt.plot(x, y, "k.")
    plt.xlabel("MJD - {:05.0f} (s)".format(mjd0))
    plt.ylabel("Frequency - {:09.0f} Hz".format(freq0))
    plt.savefig("rffit.png")
