#!/usr/bin/env python3

from astropy.time import Time

from datetime import datetime, timedelta

import argparse
import ephem
import h5py
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_artifact(filename):
    hdf5_file = h5py.File(filename, 'r')
    if hdf5_file.attrs['artifact_version'] != 2:
        print("unsupported artifact version {}".format(hdf5_file.attrs['artifact_version']))
        # return
    wf = hdf5_file.get('waterfall')
    data = (np.array(wf['data']) * np.array(wf['scale']) + np.array(wf['offset']))
    metadata = json.loads(hdf5_file.attrs['metadata'])

    return hdf5_file, wf, data, metadata

def extract_peaks(data):
    """
    Returns
    -------
    points: (time_index, freq_index)
    measurements: measurement tuples of relative time in seconds and relative frequency in hertz
    """

    snr = (np.max(data, axis=1) - np.mean(data, axis=1)) / np.std(data, axis=1)

    # Integrate each channel along the whole observation,
    # This filter helps under the assumption of minimal Doppler deviaion
    snr_integrated = np.zeros(data[0].shape)
    for row in data:
        snr_integrated += (row - np.mean(row)) / np.std(row)
    snr_integrated = pd.Series(snr_integrated/len(snr_integrated))

    SNR_CUTOFF = 4
    CHANNEL_WINDOW_SIZE = 4
    channel_mask = (snr_integrated.diff().abs() / snr_integrated.diff().std() > SNR_CUTOFF).rolling(CHANNEL_WINDOW_SIZE, min_periods=1).max()
    rows=[]
    for i, row in enumerate(channel_mask):
        if not row:
            continue
        rows.append(i)

    points_all = []
    points_filtered = []
    for i,x in enumerate(np.argmax(data, axis=1)):
        points_all.append([i, x])
        if x in rows:
            points_filtered.append([i, x])

    points_all = np.array(points_all)
    points_filtered = np.array(points_filtered)

    if len(points_filtered) == 0:
        print("WARNING: No measurement passed filter, loosening requirements now.")
        points_filtered = points_all
        measurements = None

    measurements = np.vstack((wf['relative_time'][points_filtered[:,0]],
                             [wf['frequency'][x] for x in points_filtered[:,1]]))
    return snr, measurements, points_filtered

def plot_measurements_all(wf, data):
    plt.plot(wf['relative_time'][:],
             [wf['frequency'][x] for x in np.argmax(data[:], axis=1)],
             '.',
             label='all')

def plot_measurements_selected(measurements):
    plt.plot(measurements[0,:],
             measurements[1,:],
             '.',
             label='filtered')

def plot_legend(observation_id):
    plt.legend()
    plt.title("Observation #{} - Maximum Values".format(observation_id))
    plt.grid()
    plt.xlabel('Elapsed Time / s')
    plt.ylabel('rel. Frequency / kHz')
    plt.show()



def dedoppler(measurements, metadata):
    tle = metadata['tle'].split('\n')
    start_time = datetime.strptime(wf.attrs['start_time'].decode('ascii'),
                                   '%Y-%m-%dT%H:%M:%S.%fZ')
    f_center = float(metadata['frequency'])

    # Initialize SGP4 propagator / pyephem
    satellite = ephem.readtle('sat', tle[1], tle[2])
    observer = ephem.Observer()
    observer.lat = str(metadata['location']['latitude'])
    observer.lon = str(metadata['location']['longitude'])
    observer.elevation = metadata['location']['altitude']

    def remove_doppler_correction(t, freq):
        """
        Arguments
        ---------
        t - float: Time in seconds
        freq - float: Relative Frequency in herz
        """
        observer.date = t
        satellite.compute(observer)
        v = satellite.range_velocity
        df = f_center * v / ephem.c
        return  f_center + freq - df*2

    output = []
    for dt,df in measurements.T:
        t = start_time + timedelta(seconds=dt)
        f = f_center + df
        freq_recv = remove_doppler_correction(t, f)
        output.append((t, freq_recv))
    return output

def save_rffit_data(filename, measurements, site_id):
    with open(filename, 'w') as file_out:
        for time, freq in measurements:
            line = '{:.6f}\t{:.2f}\t1.0\t{}\n'.format(Time(time).mjd, freq, site_id)
            file_out.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze SatNOGS Artifact.')
    parser.add_argument('--observation_id', type=str,
                        help='SatNOGS Observation ID')
    parser.add_argument('--filename', type=str,
                        help='SatNOGS Artifacts File')
    parser.add_argument('--output_file', type=str,
                        help='STRF-compatible output file')
    parser.add_argument('--site_id', type=int, required=True,
                        help='STRF Site ID')
    args = parser.parse_args()

    if args.observation_id:
        # Use canonic file paths, ignoring `--filename` and `--output_path`
        filename = '{}/{}.h5'.format(os.getenv('SATNOGS_ARTIFACTS_DIR'), args.observation_id)
        output_file = '{}/{}.dat'.format(os.getenv('SATNOGS_DOPPLER_OBS_DIR'), args.observation_id)
    else:
        if not any([args.filename, args.output_file]):
            print('ERROR: Missing arguments')
        filename = args.filename
        output_file = args.output_file

    print('Load {}'.format(filename))
    hdf5_file, wf, data, metadata = load_artifact(filename)

    print('Extract measurements...')
    snr, measurements, points = extract_peaks(data)
    plot_measurements_all(wf, data)
    plot_measurements_selected(measurements)
    plot_legend(args.observation_id)

    m2 = dedoppler(measurements, metadata)
    save_rffit_data(output_file, m2, site_id=args.site_id)
    print('Data written in {}'.format(output_file))
