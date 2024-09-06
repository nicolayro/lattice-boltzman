PROGRAM:=d2q6
CCLOCAL:=OMPI_CC=gcc-14 mpicc
CCPROD:=mpiicx

CFLAGS+= -std=c11 -fopenmp -Wall -Wextra -pedantic -Werror
LDLIBS+= -lm

SRC=$(PROGRAM).c ppm.c

.PHONY: build prod run images anim dirs clean purge

build: $(SRC)
	$(CCLOCAL) $^ $(CFLAGS) $(LDLIBS) -o $(PROGRAM)

prod: $(SRC)
	$(CCPROD)  $^ $(CFLAGS) $(LDLIBS) -o $(PROGRAM)

run: build
	mpirun -np 2 $(PROGRAM) -H 800 -W 1200

short: build
	mpirun -np 4 $(PROGRAM) -I 10001 -H 800 -W 1200 -F 1

geometry: build
	mpirun -np 4 $(PROGRAM) -I 10001 -s 10 -W 160 -H 40 -G name.ppm

images:
	ls data | parallel -v A={.} ./plot.sh

anim: images
	ffmpeg -y -an -i imgs/%5d.png -vcodec libx264 -pix_fmt yuv420p -profile:v baseline -level 3 -r 24 vortex_shedding.mp4

open: anim
	open vortex_shedding.mp4

dirs:
	mkdir -p data
	mkdir -p imgs

clean:
	-rm data/*.dat imgs/*.png

purge: clean
	-rm -r d2q6 data imgs
