#!/usr/bin/env python3

import argparse
import os

from satnogs_api_client import fetch_observation_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download TLE from a specific Observation in SatNOGS Network.')
    parser.add_argument('observation_id', type=int,
                        help='SatNOGS Observation ID')
    args = parser.parse_args()

    obs = fetch_observation_data([args.observation_id])[0]

    filename = '{}/{}.txt'.format(os.getenv('SATNOGS_TLE_DIR'), args.observation_id)

    with open(filename, 'w') as f:
        f.write(obs['tle0'])
        f.write('\n')
        f.write(obs['tle1'])
        f.write('\n')
        f.write(obs['tle2'])
        f.write('\n')
    print("TLE saved in {}".format(filename))
