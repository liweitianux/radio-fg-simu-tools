#
# Makefile for gen_2d tools
#
# Aaron LI
# 2015/04/14
# 2017-06-17
#

CXX= 	g++
CXXFLAGS= 	-Wall

FIO_INC= 	-I../fio/include
FIO_LIBS= 	-L../fio/lib -lfio -lcfitsio -lblitz
FIO_LDFLAGS=	-Wl,-rpath='$$ORIGIN/../fio/lib'


TARGETS= 	gen_2d_img2048 gen_2d_img_cut gen_2d_img_cut2 calc_total_flux

all: $(TARGETS)

gen_2d_img2048: gen_2d_img2048.o
	$(CXX) $(CXXFLAGS) $(FIO_LDFLAGS) $(FIO_LIBS) -o $@ $<

gen_2d_img_cut: gen_2d_img_cut.o
	$(CXX) $(CXXFLAGS) $(FIO_LDFLAGS) $(FIO_LIBS) -o $@ $<

gen_2d_img_cut2: gen_2d_img_cut2.o
	$(CXX) $(CXXFLAGS) $(FIO_LDFLAGS) $(FIO_LIBS) -o $@ $<

calc_total_flux: calc_total_flux.o
	$(CXX) $(CXXFLAGS) $(FIO_LDFLAGS) $(FIO_LIBS) -o $@ $<

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(FIO_INC) -c $<


clean:
	-rm -f *.o *~ $(TARGETS)

