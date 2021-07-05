# Guide for RF doppler analysis from SatNOGS Observation Waterfall Images

Forum post: https://community.libre.space/t/new-software-satnogs-waterfall-tabulation-helper/4380

This tool evolved from the need to tabulate doppler data from the lunar
Change-4 probe,
see <https://gitlab.com/kerel-fs/jupyter-notebooks/-/tree/master/change4/data#data-tabulation-method-new>.

## Installation
Install dependencies
```
pip install -r contrib/requirements.txt
```

# First Setup
Choose a folder where TLEs are stored and
a folder where RF doppler observatons (`.dat`-files) are stored.

For now the paths are configured via environment variables,
so make sure to set them correctly before each usage.

Example:
```
# Filename convention: {TLE_DIR}/{observation_id}.txt
SATNOGS_TLE_DIR="./data/tles"

# absulute frequency measurement storage
# Filename convention: {OBS_DIR}/{observation_id}.txt
SATNOGS_OBS_DIR="./data/obs"

# SATTOOLS/STRF/STVID sites.txt file
SATNOGS_SITES_TXT="./data/sites.txt"

mkdir -p $SATNOGS_TLE_DIR 
mkdir -p $SATNOGS_OBS_DIR 
```

## Usage
0. Make sure the (3) env variables are set.
1. Choose SatNOGS Observation ID from network and run the tabulation helper

   ```
   ./contrib/satnogs_waterfall_tabulation_helper.py 1102230
   ```

   An interactive plot will show up.
   Clicking inside the plot will add a signal marker.
   If you are finished with adding signal markers,
   save the signal markers using the keyboard shortcut `f`.
   
   Custom keyboard shortcuts:
   
   -     u - undo last signal marker
   -     f - save the signal markers in an strf-compatible file
   
   Useful Matplotlib navigation keyboard shortcuts (documentation):
   
       p - toggle 'Pan/Zoom' modus
       o - toggle 'Zoom-to-rect' modus
       h - Home/Reset (view)
       c - Back (view)

2. Run rffit for orbit fitting, e.g.
   ```
   ./rffit -d $SATNOGS_OBS_DIR/1102230.dat -i 44356 -c $SATNOGS_TLE_DIR/1102230.txt -s 7669
   ```

## Known issues
- A site id of the form `7{station_id}` is automatically assigned and written to
  the `sites.txt` (e.g. station 669 should get `7669` assigned).
  Only SatNOGS stations <999 are supported, as the strf sites.txt parse only allows
  4-digit site ids. In case of problems, choose a free site id and manually correct the
  doppler obs (`.dat`-files, last colum).
