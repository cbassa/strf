#!/usr/bin/env python3

import sys
import logging
import requests
import settings

from urllib.parse import urljoin
from pprint import pprint

from satnogs_api_client import fetch_observation_data, fetch_tle_of_observation

logger = logging.getLogger(__name__)


def fetch_latest_artifacts_metadata():
    url = urljoin(settings.SATNOGS_DB_API_URL, 'artifacts/',)
    params = {}
    headers = {'Authorization': 'Token {0}'.format(settings.SATNOGS_DB_API_TOKEN)}

    try:
        response = requests.get(url,
                                params=params,
                                headers=headers,
                                timeout=10)
        response.raise_for_status()
    except (requests.ConnectionError, requests.Timeout, requests.TooManyRedirects):
        logger.exception('An error occurred trying to GET artifact metadata from db')

    artifacts_metadata = response.json()

    if not len(artifacts_metadata):
        logger.info('No artifacts found in db')
        sys.exit(-1)

    return artifacts_metadata

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Fetch list of artifacts
    artifacts_metadata = fetch_latest_artifacts_metadata()

    observation_ids = [artifact['network_obs_id'] for artifact in artifacts_metadata]

    # Load corresponding obs from network
    observations = fetch_observation_data(sorted(observation_ids, reverse=True))

    # Filter by good status
    for observation in observations:
        if not observation['vetted_status'] == 'good':
            pass

        print("{}/observations/{}/".format(settings.SATNOGS_NETWORK_API_URL[:-5], observation['id']))
