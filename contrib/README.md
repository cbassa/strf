# STRF contrib Scripts

## Installation

0. Clone repository
   ```
   git clone https://github.com/cbassa/strf.git
   cd strf
   ```

1. Create Pyhton virtual environment
   ```
   mkvirtualenv strf
   ```

2. Install dependencies
  ```
  pip insatll -r contrib/requirements.txt
  ```

3. Create configuration based on `env-dist`
  ```
  cp env-dist .env
  ```

4. Add your SatNOGS Network API Token in the `.env` file

5. Create canonical data paths
  ```
  source .env
  mkdir -p SATNOGS_TLE_DIR
  mkdir -p SATNOGS_ARTIFACTS_DIR
  mkdir -p SATNOGS_DOPPLER_OBS_DIR
  ```

## Canonical Data Paths

For each observation there are multiple files created in different canonical locations in the
filesystem, defined by via environment variables.

| Path | Content | Example Filename |
|-|-|-|
| `SATNOGS_TLE_DIR`         | Text file containing the TLE used during an observation | `4950356.txt` |
| `SATNOGS_ARTIFACTS_DIR`   | The SatNOGS Artifacts file of an observation | `4950356.h5` |
| `SATNOGS_DOPPLER_OBS_DIR` | The STRF-compatible Frequency Measurements file | `4950356.dat` |

## Usage

0. Set environment variables
   ```
   source .env

   OBSERVATION_ID=4950356
   NORAD_ID=27844
   SITE_ID=977
   ```

1. Download TLE
   ```
   ./contrib/download_satnogs_tle.py $OBSERVATION_ID
   ```

2. Download Artifact
   ```
   ./contrib/download_satnogs_artifact.py $OBSERVATION_ID
   ```

3. Analyze Artifact
   ```
   ./contrib/analyze_artifact.py --observation_id $OBSERVATION_ID --site_id 977
   ```

4. Manually add the SatNOGS station to your site.txt

5. Fit TLE
   ```
   rffit -d "$SATNOGS_DOPPLER_OBS_DIR/$OBSERVATION_ID.dat" -c "$SATNOGS_TLE_DIR/$OBSERVATION_ID.txt" -i $NORAD_ID -s SITE_ID
   ```

In summary here the commands from above in a single listing:

```bash
OBSERVATION_ID=4950356
NORAD_ID=27844
SITE_ID=977
./contrib/download_satnogs_tle.py $OBSERVATION_ID
./contrib/download_satnogs_artifact.py $OBSERVATION_ID
./contrib/analyze_artifact.py --observation_id $OBSERVATION_ID --site_id 977
rffit -d "$SATNOGS_DOPPLER_OBS_DIR/$OBSERVATION_ID.dat" -c "$SATNOGS_TLE_DIR/$OBSERVATION_ID.txt" -i $NORAD_ID -s SITE_ID
```

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
