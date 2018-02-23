# Compiling flags
CFLAGS = -O3

# Linking flags
LFLAGS = -lcpgplot -lpgplot -lX11 -lpng -lm -lgsl -lgslcblas

# Compiler
CC = gcc

all:
	make rfedit rfplot rffft rfpng rffit

rffit: rffit.o sgdp4.o satutl.o deep.o ferror.o dsmin.o simplex.o versafit.o
	gfortran -o rffit rffit.o sgdp4.o satutl.o deep.o ferror.o dsmin.o simplex.o versafit.o $(LFLAGS)

rfpng: rfpng.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o
	gfortran -o rfpng rfpng.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

rfedit: rfedit.o rfio.o rftime.o
	$(CC) -o rfedit rfedit.o rfio.o rftime.o -lm

rfplot: rfplot.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o 
	gfortran -o rfplot rfplot.o rftime.o rfio.o rftrace.o sgdp4.o satutl.o deep.o ferror.o $(LFLAGS)

rffft: rffft.o rftime.o
	$(CC) -o rffft rffft.o rftime.o -lfftw3f -lm

clean:
	rm -f *.o
	rm -f *~
