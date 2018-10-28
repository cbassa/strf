# STRF

**strf** is the satellite tracking toolkit for radio observations (RF). The software is designed to allow tracking of satellites from radio observations, using Doppler curves to identify satellites and/or determine their orbits.

The software is designed for *linux* operating systems, and will work with most software defined radios (SDRs), certainly those that are supported by http://www.gnuradio.org. The software comes with tools for data acquisition, performing FFTs to generate timestamped spectrograms (waterfall plots), and analysis, to extract and analyse Doppler curves.

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
  * pgplot plotting library: http://www.astro.caltech.edu/~tjp/pgplot/
  * gnu scientific library: ftp://ftp.gnu.org/gnu/gsl/
  * fftw: http://www.fftw.org/download.html
* Run `make` on the **strf** folder

Run notes
---------
* You will need to set the following environment variables in your login file to run **strf**.
	* `ST_COSPAR` COSPAR number
	* `ST_DATADIR` path to **strf** directory
	* `ST_TLEDIR` path to TLE directory
	* `ST_LOGIN` space-track.org login info (of the form `ST_LOGIN="identity=username&password=password"`)
* You should install NTP support on the system and configure time/date to automatically
  synchronize to time servers.

Operation
---------
The main use of **strf** is to acquire IQ data from SDRs and produce time stamped spectrograms with the `rffft` application. `rffft` will perform Fast Fourier Transforms on the input data to a user defined number of spectral channels (via the `-c` command line option), and integrate/average these to a user defined integration length (via the `-t` command line option). The output will be a `*.bin` file which contains a 256 byte human readable header (which can be inspected with `head -c256`), followed by a binary array of floating point numbers representing the power in the spectral channels. This is an example of the 256 byte header:

	HEADER
	UTC_START    2018-01-12T15:59:13.524
	FREQ         2244000000.000000 Hz
	BW           4000000.000000 Hz
	LENGTH       0.998922 s
	NCHAN        40000
	NSUB         60
	END

The header keywords are mostly self explanatory, though the `NSUB` keyword specifies that this single `bin` file contains 60 spectra.

`rffft` can read from a previously recorded IQ recording, but is usually operated in realtime mode by reading IQ data from a so-called named pipe or fifo (first in, first out). Here, the SDR writes IQ data to a fifo (instead of a file), and `rffft` reads the samples from the fifo. Using an **airspy** as an example, it could be configured as follows:

	mkfifo fifo
	rffft -i fifo -f 101e6 -s 2.5e6 &
	airspy_rx -a 1 -f 101 -t 2 -r fifo
	
Here, we first make the fifo `mkfifo fifo`, then start `rffft` to read from the fifo (`-i` option), with a 101MHz center frequency (`-f` option) and a 2.5MHz sample rate (`-s` option). The `&` puts this command in the background. Finally, we start obtaining IQ data from the **airspy** with `airspy_rx` in the 2.5MHz sampling mode (`-a 1`) at the same frequency (`-f 101`, in MHz), with the 2.5MHz sample rate (`-t 2`) and writing the samples to the fifo (`-r fifo`). Similar scripts can be made with other SDRs, and otherwise with **gnuradio** flow graphs where the output file sink is a fifo.

The output spectrograms can be viewed and analysed using `rfplot`. 
