# Compiling flags
CFLAGS = -O3

# Linking flags
LFLAGS = -lcpgplot -lpgplot -lX11 -lpng -lm 

# Compiler
CC = gcc

all:
	make rfedit rfplot rffft

rfedit: rfedit.o
	$(CC) -o rfedit rfedit.o -lm

rfplot: rfplot.o
	gfortran -o rfplot rfplot.o $(LFLAGS)

rffft: rffft.o
	$(CC) -o rffft rffft.o -lm -lfftw3f

clean:
	rm -f *.o
	rm -f *~
