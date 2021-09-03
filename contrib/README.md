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
