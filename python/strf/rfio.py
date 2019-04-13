#!/usr/bin/env python
import numpy as np
from datetime import datetime
import re
import os
from astropy.time import Time
import scipy.optimize

class Spectrogram:
    """Spectrogram class"""

    def __init__(self, z, mjd, freq, siteid):
        """Define a spectrogram"""

        self.z = z
        self.mjd = mjd
        self.freq = freq
        self.siteid = siteid
        self.nchan = z.shape[0]
        self.nsub = z.shape[1]

    def read_bin_files(path, prefix, ifile, nsub, siteid):
        # Read first file to get number of channels
        fname = os.path.join(path, "{:s}_{:06d}.bin".format(prefix, ifile))
        with open(fname, "rb") as fp:
            header = parse_header(fp.read(256))

        # Set frequencies
        freq = np.linspace(-0.5*header["bw"], 0.5*header["bw"], header["nchan"], endpoint=False)+header["freq"]+0.5*header["bw"]/header["nchan"]
        
        # Loop over subints and files
        zs = []
        mjds = []
        isub = 0;
        while isub<nsub:
            # File name of file
            fname = os.path.join(path, "{:s}_{:06d}.bin".format(prefix, ifile))
#            print(fname)
            with open(fname, "rb") as fp:
                next_header = fp.read(256)
                while next_header:
                    header = parse_header(next_header)
                    mjds.append(Time(header["utc_start"], format="datetime", scale="utc").mjd+0.5*header["length"]/86400.0)
                    zs.append(np.fromfile(fp, dtype=np.float32, count=header["nchan"]))
                    next_header = fp.read(256)
                    isub += 1
            ifile += 1

            
        return Spectrogram(np.transpose(np.vstack(zs)), np.array(mjds), freq, siteid)

    def correct_spectrogram(self, fcen, bw):
        c = np.abs(self.freq-fcen)<bw
        x = self.freq[c]-fcen
        y = np.median(self.z, axis=1)[c]

        p = [x[np.argmax(y)], 5.0, 1.0, 0.0]

        q, cov_q = scipy.optimize.leastsq(gaussian_residual, p, args=(x, y))

        self.freq -= q[0]
    
def parse_header(header_b):
    header_s = header_b.decode('ASCII').strip('\x00')
    regex = r"^HEADER\nUTC_START    (.*)\nFREQ         (.*) Hz\nBW           (.*) Hz\nLENGTH       (.*) s\nNCHAN        (.*)\nNSUB         (.*)\nEND\n$"
    match = re.fullmatch(regex, header_s, re.MULTILINE)

    utc_start = datetime.strptime(match.group(1), '%Y-%m-%dT%H:%M:%S.%f')

    return {'utc_start': utc_start,
            'freq': float(match.group(2)),
            'bw': float(match.group(3)),
            'length': float(match.group(4)),
            'nchan': int(match.group(5)),
            'nsub': int(match.group(6))}
    
def gaussian(x, a):
    arg = -0.5*((x-a[0])/a[1])**2
    return a[2]*np.exp(arg)+a[3]

def gaussian_residual(a, x, y):
    ym = gaussian(x, a)
    return y-ym
