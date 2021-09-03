#!/usr/bin/env python3

import argparse
import tempfile
import sys
import logging
import requests
import settings

from urllib.parse import urljoin
from pprint import pprint


logger = logging.getLogger(__name__)

def fetch_artifact_metadata(network_obs_id):
    url = urljoin(settings.SATNOGS_DB_API_URL, 'artifacts/',)
    params = {'network_obs_id': network_obs_id}
    headers = {'Authorization': 'Token {0}'.format(settings.SATNOGS_DB_API_TOKEN)}

    response = requests.get(url,
                            params=params,
                            headers=headers,
                            timeout=10)
    response.raise_for_status()

    return response.json()


def fetch_artifact(url, artifact_filename):
    headers = {'Authorization': 'Token {0}'.format(settings.SATNOGS_DB_API_TOKEN)}

    response = requests.get(url,
                            headers=headers,
                            stream=True,
                            timeout=10)
    response.raise_for_status()

    with open(artifact_filename, 'wb') as fname:
        for chunk in response.iter_content(chunk_size=1024):
            fname.write(chunk)


def download_artifact(observation_id):
    try:
        artifact_metadata = fetch_artifact_metadata(network_obs_id=observation_id)
    except requests.HTTPError:
        print('An error occurred trying to GET artifact metadata from db')
        return

    if not len(artifact_metadata):
        print('No artifact found in db for network_obs_id {}'.format(network_obs_id))
        return

    print("Artifact Metadata for Observation #{} found.".format(observation_id))

    try:
        artifact_file_url = artifact_metadata[0]['artifact_file']
        artifact_file = tempfile.NamedTemporaryFile(delete=False)
        fetch_artifact(artifact_file_url, artifact_file.name)
        print("Artifact for Observation #{} saved in '{}'".format(observation_id, artifact_file.name))
    except requests.HTTPError:
        print('Download failed for {}'.format(artifact_file_url))
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download SatNOGS Artifacts from SatNOGS DB.')
    parser.add_argument('observation_ids', metavar='ID', type=int, nargs='+',
                        help='SatNOGS Observation ID')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    for observation_id in args.observation_ids:
        download_artifact(observation_id)
