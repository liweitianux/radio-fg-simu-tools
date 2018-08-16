target=rescale fill_hole

all:$(target)

rescale:rescale.cc
	g++ $< -o $@ -I /usr/local/src/opt_utilities/ -lfio -lcfitsio -I /usr/include/gsl -lgsl -lgslcblas -g -O3


fill_hole:fill_hole.cc
	g++ $< -o $@ -I /usr/local/src/opt_utilities/ -lfio -lcfitsio -I /usr/include/gsl -lgsl -lgslcblas -g -O3


clean:
	rm -f *.o $(target)
