# STRF contrib Scripts

## Installation

- Create Pyhton virtual environment and install dependencies via pip:
  ```
  mkvirtualenv strf
  pip install -r contrib/requirements.txt
  ```

- Create initial configuration from `env-dist`
  ```
  cp env-dist .env
  ```
- Configure the data paths (they must already exist!) and add your SatNOGS DB API Token in `.env`

## Find good observations with SatNOGS artifacts
```
$ ./contrib/find_good_satnogs_artifacts.py
```

## Download SatNOGS Artifacts

```
$ ./contrib/download_satnogs_artifact.py 4950356
Artifact Metadata for Observation #4950356 found.
Artifact saved in /home/pi/data/artifacts/4950356.h5
```

## Download SatNOGS Observation TLEs

To add the TLE used in a specific SatNOGS Observation to your catalog, the following command can be used:
```
$ ./contrib/download_satnogs_tle.py 4950356
TLE saved in /home/pi/data/tles/satnogs/4950356.txt
```

## Plot SatNOGS Artifacts

This can be done using [spectranalysis](https://github.com/kerel-fs/spectranalysis).

## Ultra-Short Analysis Guide

```
source .env

OBSERVATION_ID=4950356
./contrib/download_satnogs_tle.py $OBSERVATION_ID
./contrib/download_satnogs_artifact.py $OBSERVATION_ID
./contrib/analyze_artifact.py --observation_id $OBSERVATION_ID --site_id 977
rffit -d "$SATNOGS_DOPPLER_OBS_DIR/$OBSERVATION_ID.dat" -c "$SATNOGS_TLE_DIR/$OBSERVATION_ID.txt" -i 27844 -s 977
```

```
$ ./contrib/download_satnogs_tle.py 4950356
TLE saved in /mnt/old_home/kerel/c4/satnogs/data/tles/satnogs/4950356.txt
$ ./contrib/download_satnogs_artifact.py 4950356
Artifact Metadata for Observation #4950356 found.
Artifact saved in /mnt/old_home/kerel/c4/satnogs/data/artifacts/4950356.h5
$ ./contrib/analyze_artifact.py --observation_id 4950356 --site_id 977
Load /mnt/old_home/kerel/c4/satnogs/data/artifacts/4950356.h5
Extract measurements...
Data written in /mnt/old_home/kerel/c4/satnogs/data/doppler_obs/4950356.dat
$ rffit -d ${SATNOGS_DOPPLER_OBS_DIR}/4950356.dat -c ${SATNOGS_TLE_DIR}/4950356.txt -i 27844 -s 977
```
