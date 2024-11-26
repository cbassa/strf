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

rffit: rffit.o sgdp4.o satutl.o deep.o ferror.o dsmin.o simplex.o versafit.o rfsites.o
	gfortran -o rffit rffit.o sgdp4.o satutl.o deep.o ferror.o dsmin.o simplex.o versafit.o rfsites.o $(LFLAGS)

rfpng: rfpng.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o rftles.o
	gfortran -o rfpng rfpng.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o rftles.o $(LFLAGS)

rfdop: rfdop.o rftrace.o rfio.o rftime.o sgdp4.o satutl.o deep.o ferror.o rftles.o
	$(CC) -o rfdop rfdop.o rftrace.o rfio.o rftime.o sgdp4.o satutl.o deep.o ferror.o rftles.o -lm

rfedit: rfedit.o rfio.o rftime.o
	$(CC) -o rfedit rfedit.o rfio.o rftime.o -lm

rffind: rffind.o rfio.o rftime.o
	$(CC) -o rffind rffind.o rfio.o rftime.o -lm

rftrack: rftrack.o rfio.o rftime.o rftrace.o sgdp4.o satutl.o deep.o ferror.o
	$(CC) -o rftrack rftrack.o rfio.o rftime.o rftrace.o sgdp4.o satutl.o deep.o ferror.o -lm

rfplot: rfplot.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o rftles.o
	gfortran -o rfplot rfplot.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o versafit.o dsmin.o simplex.o rftles.o $(LFLAGS)

rffft: rffft.o rffft_internal.o rftime.o
	$(CC) -o rffft rffft.o rffft_internal.o rftime.o -lfftw3f -lm -lsox

tests/tests: tests/tests.o tests/tests_rffft_internal.o tests/tests_rftles.o rffft_internal.o rftles.o satutl.o ferror.o
	$(CC) -Wall -o $@ $^ -lcmocka -lm

tests: tests/tests
	./tests/tests

.PHONY: clean install uninstall tests

clean:
	rm -f *.o tests/*.o
	rm -f *~

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
