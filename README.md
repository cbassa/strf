# STRF

**strf** is the satellite tracking toolkit for radio observations (RF). The software is designed to allow tracking of satellites from radio observations, using Doppler curves to identify satellites and/or determine their orbits.

The software is designed for *linux* operating systems, and will work with most software defined radios (SDRs), certainly those that are supported by http://www.gnuradio.org. The software comes with tools for data acquisition, performing FFTs to generate timestamped spectrograms (waterfall plots), and analysis, to extract and analyse Doppler curves.

Install
------

* For Ubuntu systems or similar.
  * Install dependencies: `sudo apt install git make gcc pgplot5 gfortran libpng-dev libx11-dev libgsl-dev libfftw3-dev dos2unix`
  * Clone repository: `git clone https://github.com/cbassa/strf.git`
  * Compile: `cd strf; make`
  * Install (in `/usr/local`): `sudo make install`

Configure
---------
* You will need to set the following environment variables in your login file to run **strf**.
	* `ST_DATADIR` path to **strf** directory (e.g. `$HOME/software/strf`)
	* `ST_TLEDIR` path to TLE directory (e.g. `$HOME/tle`)
	* `ST_COSPAR` COSPAR site number (add to site location to `$ST_DATADIR/data/sites.txt`)
	* `ST_LOGIN` space-track.org login info (of the form `ST_LOGIN="identity=username&password=password"`)
* Run `tleupdate` to download latest TLEs.
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

Alternatively, when no input filename is given (with the `-i` option), `rffft` will read from stdin so it is possible to directly pipe an SDR receiver's application into `rffft`.

With an RTL-SDR:

    rtl_sdr -g 29 -f 97400000 -s 2048000 - | ./rffft -f 97400000 -s 2048000 -F char

Here we use the **RTL-SDR** receiver with `rtl_sdr` with a gain of 29dB (`-g 29`), a center frequency of 97.4MHz (`-f 97400000`, in Hz) and a samplerate of 2.048MS/s (`-s 2048000`, in S/s). Note the trailing dash (`-`) in the `rtl_sdr` command to tell it to write to stdout instead of a file so it can be piped (`|`) through `rffft`. The same center frequency and samplerate are given to `rffft`. As `rtl_sdr` outputs data as 8 bits, `-F char` is required to tell `rffft` the format of the data.

With a HackRF:

    hackrf_transfer -l 24 -g 32 -f 97400000 -s 8000000 -r - | ./rffft -f 97400000 -s 8000000 -F char -c 100

Here we use the **HackRF** receiver with `hackrf_transfer` with a lna gain of 24dB (`-l 24`), an IF gain of 32dB (`-g 32`), a center frequency of 97.4MHz (`-f 97400000`, in Hz) and a samplerate of 8MS/s (`-s 8000000`, in S/s). The output file is given as stdout (`-r -`). Again the same frequency and samplerate are given to `rffft` and as `hackrf_transfer` also outputs 8 bit data `-F char` is also required for `rffft`.

With a Adalm Pluto:

    iio_attr -u usb:x.y.z -c ad9361-phy RX_LO frequency 97400000
    iio_attr -u usb:x.y.z -c ad9361-phy voltage0 rf_port_select A_BALANCED
    iio_attr -u usb:x.y.z -c ad9361-phy voltage0 rf_bandwidth 2000000
    iio_attr -u usb:x.y.z -c ad9361-phy voltage0 sampling_frequency 2000000
    iio_attr -u usb:x.y.z -c ad9361-phy voltage0 gain_control_mode manual
    iio_attr -u usb:x.y.z -c ad9361-phy voltage0 hardwaregain 60
    iio_readdev -u usb:x.y.z -b 4194304 cf-ad9361-lpc | ./rffft -f 97400000 -s 2000000 -F int

Here we use the **Adalm Pluto** transceiver with `iio_readdev`. The transceiver is connected via USB. Replace `x.y.z` by the USB bus ID of your transceiver. You can retreive the USB IDs from the output of `iio_info -s`. Connection via Ethernet is possible as well. Before receiving we have to set center frequency, samplerate, gain, ... by a bunch of `iio_attr` commands.

With I/Q recordings obtained from Gqrx:

    ./rffft -i gqrx_YYYYMMDD_HHMMSS_97400000_2000000_fc.raw -f 97400000 -s 2000000 -F float -T "YYYY-MM-DDTHH:MM:SS"

**Gqrx** records complex samples into `raw` files. The filename contains date, time, center frequency and samplerate separated by underscores. Replace `YYYYMMDD` and `HHMMSS` by your actual time and respectively. Pay attention to insert an uppercase `T` between date and time in the time stamp parameter of the `rffft` command.

The output spectrograms can be viewed and analysed using `rfplot`. 
