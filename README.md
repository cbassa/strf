# STRF

**strf** is the satellite tracking toolkit for radio observations (RF). The software is designed to allow tracking of satellites from radio observations, using Doppler curves to identify satellites and/or determine their orbits.

The software is designed for *linux* operating systems, and will work with most software defined radios (SDRs), certainly those that are supported by http://www.gnuradio.com. The software comes with tools for data acquisition, performing FFTs to generate timestamped spectrograms (waterfall plots), and analysis, to extract and analyse Doppler curves.

Install
------
* Clone locally the code repository
* Install common dependencies
  * gfortran
  * gcc
  * libpng-dev
  * libx11-dev
  * libjpeg-dev
  * libexif-dev
* Build & install required libraries
  * pgplot-5.2.2: http://www.astro.caltech.edu/~tjp/pgplot/
  * gsl-2.4: ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
* Run `make` on the **strf** folder

Run notes
---------
* You will need to set the following environment variables to run **strf**.
	`ST_COSPAR` COSPAR number
	`ST_DATADIR` path to **strf** directory
	`ST_TLEDIR` path to TLE directory
	`ST_OBSDIR` path to observations directory
* You should install NTP support on the system and configure time/date to automatically
  sinchronize to time servers.
