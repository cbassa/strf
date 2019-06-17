from datetime import datetime
import numpy as np
import glob
import os
import re
from pathlib import Path


def read_spectrum(path):
    # Get the number of files
    filenames = glob.glob(os.path.join(path, '*.bin'))
    datestr = Path(filenames[0]).stem.split('_')[0]

    # Read first file to get the number of channels and number of "samples"
    filename = os.path.join(path, '{:s}_{:06}.bin'.format(datestr, 0))
    with open(filename, 'rb') as f:
        header = parse_header(f.read(256))

    zs = []
    headers = []
    for i_file in range(len(filenames)):
        filename = os.path.join(path, '{:s}_{:06}.bin'.format(datestr, i_file))
        # i_sub = 0
        with open(filename, 'rb') as f:
            next_header = f.read(256)
            while(next_header):
                headers.append(parse_header(next_header))
                zs.append(np.fromfile(f, dtype=np.float32, count=header["nchan"]))
                next_header = f.read(256)
    return np.transpose(np.vstack(zs)), headers


def parse_header(header_b):
    # TODO. Support files with the additional fields
    #       - NBITS
    #       - MEAN
    #       - RMS
    # "HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\n"

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
