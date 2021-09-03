# STRF contrib Scripts

## Installation

- Create Pyhton virtual environment and install dependencies via pip:
  ```
  mkvirtualenv strf
  pip insatll -r contrib/requirements.txt
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
$ ./contrib/download_satnogs_artifact.py 2786575
Artifact Metadata for Observation #2786575 found.
Download failed for https://db-satnogs.freetls.fastly.net/media/artifacts/b4975058-04eb-4ab7-9c40-9bcce76d94db.h5
```

```
$ ./contrib/download_satnogs_artifact.py 4443137
Artifact Metadata for Observation #4443137 found.
Artifact for Observation #4443137 saved in '/tmp/tmp1rdpnz_k'
```
