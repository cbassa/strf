# SGP4/SDP4 library (git submodule) and strf's additions to it
SGP4DIR = extern/dundee_sgp4
SGP4LIB = $(SGP4DIR)/libsgdp4.a
SGP4 = strf_util.o $(SGP4LIB)
SGP4HDRS = $(SGP4DIR)/sgdp4h.h $(SGP4DIR)/satutl.h $(SGP4DIR)/elset.h $(SGP4DIR)/alpha5.h $(SGP4DIR)/omm.h

# Library include path, kept out of CFLAGS so that overriding CFLAGS keeps it
override CPPFLAGS += -I$(SGP4DIR)

# Compiling flags
CFLAGS = -O3

# Linking flags
LFLAGS = -lcpgplot -lpgplot -lX11 -lpng -lm -lgsl -lgslcblas

# Compiler
CC = gcc

# Installation
INSTALL_PROGRAM = install -m 755
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin

all:
	make rfedit rfplot rffft rfpng rffit rffind rfdop

rffit: rffit.o dsmin.o simplex.o versafit.o rfsites.o rftles.o $(SGP4)
	gfortran -o rffit rffit.o dsmin.o simplex.o versafit.o rfsites.o rftles.o $(SGP4) $(LFLAGS)

rfpng: rfpng.o rftime.o rfio.o rftrace.o rftles.o zscale.o $(SGP4)
	gfortran -o rfpng rfpng.o rftime.o rfio.o rftrace.o rftles.o zscale.o $(SGP4) $(LFLAGS)

rfdop: rfdop.o rftrace.o rfio.o rftime.o rftles.o zscale.o $(SGP4)
	$(CC) -o rfdop rfdop.o rftrace.o rfio.o rftime.o rftles.o zscale.o $(SGP4) -lm

rfedit: rfedit.o rfio.o rftime.o zscale.o
	$(CC) -o rfedit rfedit.o rfio.o rftime.o zscale.o -lm

rffind: rffind.o rfio.o rftime.o zscale.o
	$(CC) -o rffind rffind.o rfio.o rftime.o zscale.o -lm

rftrack: rftrack.o rfio.o rftime.o rftrace.o zscale.o $(SGP4)
	$(CC) -o rftrack rftrack.o rfio.o rftime.o rftrace.o zscale.o $(SGP4) -lm

rfplot: rfplot.o rftime.o rfio.o rftrace.o versafit.o dsmin.o simplex.o rftles.o zscale.o $(SGP4)
	gfortran -o rfplot rfplot.o rftime.o rfio.o rftrace.o versafit.o dsmin.o simplex.o rftles.o zscale.o $(SGP4) $(LFLAGS)

rffft: rffft.o rffft_internal.o rftime.o
	$(CC) -o rffft rffft.o rffft_internal.o rftime.o -lfftw3f -lm -lsox

tests/tests: tests/tests.o tests/tests_rffft_internal.o tests/tests_rfio.o tests/tests_rftles.o rffft_internal.o rfio.o rftime.o rftles.o strf_util.o zscale.o $(SGP4LIB)
	$(CC) -Wall -o $@ $^ -lcmocka -lm

# SGP4/SDP4 library and element set readers, from the dundee_sgp4 submodule
$(SGP4LIB): FORCE
	$(MAKE) -C $(SGP4DIR) lib CC="$(CC)"

# Every object depends on the library headers (rebuilt when the submodule
# changes); a missing header means the submodule is not checked out
$(patsubst %.c,%.o,$(wildcard *.c tests/*.c)): $(SGP4HDRS)

$(SGP4HDRS):
	@echo "$(SGP4DIR) is empty: run 'git submodule update --init'"; exit 1

FORCE:

tests/random:
	dd if=/dev/random of=tests/random bs=12288000 count=4

tests/bins: rffft tests/random
	mkdir -p tests/bin/
	./rffft -i tests/random -p tests/bin/ -o t10_n120_32 -f 2250000000 -s 1024 -T 2016-09-28T00:20:00 -t 10 -n 120
	./rffft -i tests/random -p tests/bin/ -o t15_n90_8 -f 2250000000 -s 1024 -T 2016-09-28T00:20:00 -t 15 -n 90 -b
	touch tests/bins

tests: tests/tests tests/bins
	./tests/tests

.PHONY: clean install uninstall tests FORCE

clean:
	rm -f *.o tests/*.o
	rm -f *~
	-$(MAKE) -C $(SGP4DIR) clean

install:
	$(INSTALL_PROGRAM) rffit $(DESTDIR)$(bindir)/rffit
	$(INSTALL_PROGRAM) rfpng $(DESTDIR)$(bindir)/rfpng
	$(INSTALL_PROGRAM) rfedit $(DESTDIR)$(bindir)/rfedit
	$(INSTALL_PROGRAM) rffind $(DESTDIR)$(bindir)/rffind
	$(INSTALL_PROGRAM) rfplot $(DESTDIR)$(bindir)/rfplot
	$(INSTALL_PROGRAM) rffft $(DESTDIR)$(bindir)/rffft
	$(INSTALL_PROGRAM) tleupdate $(DESTDIR)$(bindir)/tleupdate

uninstall:
	$(RM) $(DESTDIR)$(bindir)/rffit
	$(RM) $(DESTDIR)$(bindir)/rfpng
	$(RM) $(DESTDIR)$(bindir)/rfedit
	$(RM) $(DESTDIR)$(bindir)/rffind
	$(RM) $(DESTDIR)$(bindir)/rfplot
	$(RM) $(DESTDIR)$(bindir)/rffft
	$(RM) $(DESTDIR)$(bindir)/tleupdate
