##
## libfio Makefile
## C++ bindings to libcfitsio
##
## Junhua Gu
## v1.0, 2007/06/07
##
## Weitian LI
## v1.1, 2012-05-08
## v1.2, 2017-06-17
## v1.3, 2018-01-02
##


CXX?= 		g++
CXXFLAGS= 	-g -Wall -pipe -fPIC -I./include
LDFLAGS=	-L./lib -Wl,-rpath='$$ORIGIN/lib' -lfio -lblitz -lcfitsio

OBJS= 		fio.o fits_trait.o fitsfile.o region.o region_imp.o
TGT= 		lib/libfio.so include/fio.h fits2txt txt2fits


all: $(TGT)


lib/libfio.so: libfio.o $(OBJS)
	-mkdir lib
	$(CXX) $(CXXFLAGS) -shared -o $@ $< $(OBJS)

include/fio.h: fio.h
	-mkdir include
	cp $< $@


fits2txt: fits2txt.o lib/libfio.so lib/libblitz.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<


txt2fits: txt2fits.o lib/libfio.so lib/libblitz.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<


%.o: %.cc blitz
	$(CXX) $(CXXFLAGS) -c $<


blitz: lib/libblitz.so


lib/libblitz.so: blitz-0.9
	-mkdir lib
	cd $< && \
	    ./configure --prefix=/tmp/$<-install \
	        --libdir=$(PWD)/lib --include=$(PWD)/include \
	        --disable-fortran --disable-doxygen \
	        --disable-dot --disable-html-docs \
	        --disable-static --enable-shared; \
	    make && \
	    make install && \
	    rm -rf /tmp/$<-install


blitz-0.9: blitz-0.9.tar.gz
	tar xf $<

clean:
	rm -f *.o $(TGT)
	rm -rf blitz-0.9 include/blitz include/random
	rm -rf lib/pkgconfig lib/libblitz.*


.PHONY: clean blitz

#  vim: set ts=8 sw=4 tw=0 fenc=utf-8 ft=make: #
