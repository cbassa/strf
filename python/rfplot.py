#!/usr/bin/env python
import numpy as np
from datetime import datetime
import glob
import os
import re

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


def read_spectrogram(path, prefix, isub, nsub):
    fname = os.path.join(path, "{:s}_{:06d}.bin".format(prefix, isub))
    
    # Read first file to get number of channels
    with open(fname, "rb") as fp:
        # Read header
        header = parse_header(fp.read(256))

    nchan = header["nchan"]
    print(nchan)
        
if __name__ == "__main__":
    # Settings
    path = "/data2/satobs/radio/graves_dt/bin"
    prefix = "2018-12-31T09:34:11"

    read_spectrogram(path, prefix, 48000, 900)
