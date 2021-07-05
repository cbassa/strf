#!/usr/bin/env python3

import argparse
import requests
import ephem
import csv


from astropy.time import Time
from astropy import constants as const
from datetime import datetime, timedelta
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import settings
from satnogs_api_client import fetch_observation_data, fetch_tle_of_observation


def tabulation_helper_dialog(lower_left, upper_right,
                             bandwidth, duration,
                             f_center, filename_out,
                             waterfall_matrix,
                             epoch_start, site_id, correction_method=None):
    # Assumptions: x-axis - freqency, y-axis - time; up-time advance, right-freq advance

    # Derive conversion parameters
    x_range = upper_right[0] - lower_left[0]
    y_range = lower_left[1] - upper_right[1]
    x_offset = int(lower_left[0] + 0.5 * x_range)
    y_offset = lower_left[1]

    t_step = duration / y_range
    f_step = bandwidth / x_range

    # Highlight center and calibration lines
    waterfall_matrix[:,x_offset] = (255,0,0,255)
    waterfall_matrix[y_offset,:] = (255,0,0,255)
    waterfall_matrix[:,upper_right[0]] = (255,0,0,255)
    waterfall_matrix[upper_right[1],:] = (255,0,0,255)


    markers = {'x': [], 'y': [], 't': [], 'f': []}
    
    fig, ax = plt.subplots()
    ax.imshow(waterfall_matrix)
    marker_line, = ax.plot([0], [0], marker='.', c='k', zorder=100)
    def on_mouse_press(event):
        tb = plt.get_current_fig_manager().toolbar
        if not (event.button==1 and event.inaxes and tb.mode == ''):
            return
    
        markers['x'].append(event.xdata)
        markers['y'].append(event.ydata)
        markers['t'].append(epoch_start - (event.ydata - y_offset) * t_step)
        markers['f'].append(f_step * (event.xdata - x_offset))
        print('x:{:.2f}  y:{:.2f}  --> {:7.0f} -->  {}   {:7.0f}'.format(markers['x'][-1],
                                                                         markers['y'][-1],
                                                                         markers['f'][-1] - f_center,
                                                                         markers['t'][-1],
                                                                         markers['f'][-1]))
        update_plot_markers()

    def update_plot_markers():
        marker_line.set_xdata(markers['x'])
        marker_line.set_ydata(markers['y'])
        # ax.scatter(x=[int(event.xdata)], y=[int(event.ydata)], marker='.', c='k', zorder=100)
        ax.figure.canvas.draw()

    def undo_track_point_selection():
        markers['x'] = markers['x'][:-1]
        markers['y'] = markers['y'][:-1]
        markers['t'] = markers['t'][:-1]
        markers['f'] = markers['f'][:-1]
        update_plot_markers()

    def write_strf_spectrum(filename, site_id, markers):
        with open(filename, 'w') as f:
            for t,freq in list(sorted(zip(markers['t'],markers['f']), key=lambda x: x[0])):

                if correction_method:
                    freq_recv = correction_method(t, freq)
                else:
                    freq_recv = freq
                line = '{:.6f}\t{:.2f}\t1.0\t{}\n'.format(Time(t).mjd, freq_recv, site_id)
                print(line, end='')
                f.write(line)

    def save_selected_track_points():
        write_strf_spectrum(filename_out, site_id, markers)
        print('Stored {} selected track points in {}.'.format(len(markers['t']), filename_out))
    
    def on_key(event):
        if event.key == 'u':
            undo_track_point_selection()
        elif event.key == 'f':
            save_selected_track_points()
    
    fig.canvas.mpl_connect('button_press_event', on_mouse_press)
    fig.canvas.mpl_connect('key_press_event', on_key)
    plt.tight_layout()
    plt.show()


def read_sitestxt(path):
    sites = {}
    with open(path) as fp:
        d = csv.DictReader([row for row in fp if not row.startswith('#')],
                           delimiter=' ',
                           fieldnames=["no", "id", "lat", "lon", "alt"],
                           restkey="observer",
                           skipinitialspace=True)
        for line in d:
            print(line)
            sites[line["no"]] = {'lat': line["lat"],
                                 'lon': line["lon"],
                                 'sn': line["id"],
                                 'alt': line["alt"],
                                 'name': ' '.join(line["observer"])}
    return sites


def add_station_to_sitestxt(station):
    # Read the sites.txt
    sites = read_sitestxt(settings.SITES_TXT)
    if '7{}'.format(station['id']) in sites:
        # Station is already available
        return
    else:
        with open(settings.SITES_TXT, 'a') as f:
            f.write('7{} {} {} {} {} {}\n'.format(station['id'],
                                                   'SN',
                                                   station['lat'],
                                                   station['lng'],
                                                   station['alt'],
                                                   station['name']))


def tabulation_helper(observation_id):
    # Fetch waterfall image and observation data
    observation = fetch_observation_data([observation_id])[0]
    tle = fetch_tle_of_observation(observation_id)
    result = requests.get(observation['waterfall'])
    image = Image.open(BytesIO(result.content))
    waterfall_matrix = np.array(image)

    epoch_start = datetime.strptime(observation['start'], '%Y-%m-%dT%H:%M:%SZ')
    epoch_end = datetime.strptime(observation['end'], '%Y-%m-%dT%H:%M:%SZ')
    site_id = int("7" + str(observation['ground_station']))

    station = {'lat': observation['station_lat'],
               'lng': observation['station_lng'],
               'alt': observation['station_alt'],
               'name': "SatNOGS No.{}".format(observation['ground_station']),
                'id': observation['ground_station']}

    # Store TLE to file
    with open('{}/{}.txt'.format(settings.TLE_DIR, observation_id), 'w') as f:
        f.write('0 OBJECT\n')
        for line in tle:
            f.write(f'{line}\n')

    # Add SatNOGS station to sites.txt
    add_station_to_sitestxt(station)

    # Set image parameters
    lower_left = (66, 1553)
    upper_right = (686, 13)

    # Set data parameters
    bandwidth = 48e3 # Hz
    duration = epoch_end - epoch_start
    f = observation['transmitter_downlink_low']
    drift = observation['transmitter_downlink_drift']
    if drift:
        f_center = f * (1 + drift / 1e9)
    else:
        f_center = f
    filename_out = '{}/{}.dat'.format(settings.OBS_DIR, observation_id)

    # Initialize SGP4 propagator / pyephem
    satellite = ephem.readtle('sat', tle[0], tle[1])
    observer = ephem.Observer()
    observer.lat = str(station['lat'])
    observer.lon = str(station['lng'])
    observer.elevation = station['alt']


    def remove_doppler_correction(t, freq):
        # Remove the doppler correction
        observer.date = t
        satellite.compute(observer)
        v = satellite.range_velocity
        df = f_center * v / const.c.value
        return  f_center + freq - df

    tabulation_helper_dialog(lower_left, upper_right,
                             bandwidth, duration,
                             f_center, filename_out,
                             waterfall_matrix,
                             epoch_start, site_id, correction_method=remove_doppler_correction)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Interactive helper to tabulate signals from SatNOGS waterfalls.')
    parser.add_argument('observation_ids', metavar='ID', type=int, nargs='+',
                        help='SatNOGS Observation ID')

    args = parser.parse_args()

    for observation_id in args.observation_ids:
        tabulation_helper(observation_id)
