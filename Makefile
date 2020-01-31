CC=gcc -O3
LDFLAGS=-lm -lgmp -fopenmp
CFLAGS=-D CONST_L=$L -D ADDRESS_NUM_BITS=$A -D BITMAP_NUM_BITS=$B
SOURCES=common.c
OBJECTS=$(SOURCES:.c=.o)

all: clean precomputation main

precomputation: precomputation.c
	$(CC) precomputation.c $(SOURCES) $(LDFLAGS) $(CFLAGS) -o precomp.out

main: main.c
	$(CC) main.c $(SOURCES) $(LDFLAGS) $(CFLAGS) -o leg.out

clean:
	rm -f *.out