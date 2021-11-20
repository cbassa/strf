# very opioninated configuration file by and for kerel

# KEEP IN SYNC WITH stvid/configuration.ini !!!!

# export PGPLOT_DIR=/usr/local/src/pgplot-5.2.2
export PGPLOT_FONT=/usr/lib/grfont.dat

source ../strf/.env

SATTOOLS=$SATNOGS_DIR/sattools

export ST_COSPAR
export ST_DATADIR=$SATNOGS_DIR/data/
export ST_TLEDIR=$SATNOGS_DIR/data/tles/sattools
export ST_OBSDIR=$SATNOGS_DIR/data/satobs

# Only on dev branch:
export ST_SITES=$SATNOGS_DIR/data/data/sites.txt

# space-track.org credentials
export ST_LOGIN

PATH=$PATH:$SATTOOLS
export PATH=$PATH:$SATTOOLS/scripts
