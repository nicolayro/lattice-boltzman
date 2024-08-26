CFLAGS+= -std=c99 -O2 -Wall -Wextra -pedantic -Werror
LDLIBS+= -lm
TARGETS= lbm_vortex dirs
IMAGES=$(shell ls data/*.dat | sed s/data/imgs/g | sed s/\.dat/.png/g)
.PHONY: all images anim dirs clean purge
all: ${TARGETS}
images: ${IMAGES}
anim: images
	ffmpeg -y -an -i imgs/%5d.png -vcodec libx264 -pix_fmt yuv420p -profile:v baseline -level 3 -r 12 vortex_shedding.mp4
dirs:
	mkdir -p data
	mkdir -p imgs
imgs/%.png: data/%.dat
	echo "set term png size 800,600; set output \"imgs/$*.png\"; set view 0,0,1; set cbrange [0:0.6]; set palette defined (0 \"black\",12 \"cyan\", 16\"white\"); splot \"data/$*.dat\" binary matrix with pm3d" | gnuplot -
clean:
	-rm *.o data/*.dat imgs/*.png
purge: clean
	-rm -r lbm_vortex data imgs
