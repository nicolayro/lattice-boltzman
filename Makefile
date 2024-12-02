PROGRAM:=d2q6
CCLOCAL:=OMPI_CC=gcc-14 mpicc
CCPROD:=mpiicx

CFLAGS+= -std=c11 -fopenmp -Wall -Wextra -pedantic -Werror -O2
LDLIBS+= -lm

SRC=$(PROGRAM).c ppm.c

THREADS=1
PROCS=1

.PHONY: build prod run images anim dirs clean purge

build: $(SRC)
	$(CCLOCAL) $^ $(CFLAGS) $(LDLIBS) -o $(PROGRAM)

prod: $(SRC)
	$(CCPROD)  $^ $(CFLAGS) $(LDLIBS) -o $(PROGRAM)

run: build
	OMP_NUM_THREADS=$(THREADS) mpirun -np $(PROCS) $(PROGRAM) -I 1000 -s 1 -G assets/circle.ppm

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

7: d2q7.c
	$(CCLOCAL) $^ $(CFLAGS) $(LDLIBS) -o d2q7

run_d2q7: 7
	mpirun -np $(PROCS) d2q7 -i 40000 assets/circle.ppm
