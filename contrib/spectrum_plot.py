from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def plot_spectrum(z, headers):
    z_mean = np.mean(z)
    z_std = np.std(z)
    zmin = z_mean - 2 * z_std
    zmax = z_mean + 6 * z_std
    
    f_min = -headers[0]['bw']*1e-3
    f_max = +headers[0]['bw']*1e-3
    
    utc_start = headers[0]['utc_start']
    t_min = mdates.date2num(utc_start + timedelta(seconds=0))
    t_max = mdates.date2num(utc_start + timedelta(seconds=z.shape[1]))
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.imshow(z,
              origin='lower',
              aspect='auto',
              extent=[t_min, t_max, f_min, f_max],
              vmin=zmin,
              vmax=zmax,
              interpolation='None',
              cmap='viridis')
    
    ax.xaxis_date()
    date_format = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(date_format)
    
    # diagonal xaxis labels
    fig.autofmt_xdate()
    
    plt.xlabel('Time;  start: {:%Y-%m-%d %H:%M:%S}Z'.format(headers[0]['utc_start']))
    plt.ylabel('Freqeuncy / kHz; center: {} MHz'.format(headers[0]['freq'] * 1e-6))

    plt.tight_layout()

    plt.show()
