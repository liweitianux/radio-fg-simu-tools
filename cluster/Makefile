#
# Makefile for the tools used to simulate and generate cluster images.
#
# Aaron LI
# 2015/04/12
#

CXX= 	g++
CXXFLAGS= 	-Wall

CC= 	gcc
CFLAGS= 	-Wall

GSL_LIBS= 	-lgsl -lgslcblas -lm

FIO_INC=	-I./fio/include
FIO_LIBS= 	-L./fio/lib -lfio -lcfitsio
FIO_LDFLAGS=	-Wl,-rpath='$$ORIGIN/fio/lib'

TARGETS= 	sphxy alpha_sim beta_sim beta_calc calc_totalnum \
		calc_T_P1400 calc_s0 calc_F1400_rcDA \
		z_m_simulation transfer_sur \
		calc_nmz draw_line histogram \
		calc_NN0 calc_Fnu_t

all: $(TARGETS)

sphxy: sphxy.o
	$(CXX) $(CXXFLAGS) -o $@ $<

alpha_sim: alpha_sim.o
	$(CXX) $(CXXFLAGS) -o $@ $<

beta_sim: beta_sim.o
	$(CXX) $(CXXFLAGS) -o $@ $<

beta_calc: beta_calc.o
	$(CXX) $(CXXFLAGS) -o $@ $<

calc_totalnum: calc_totalnum.o
	$(CXX) $(CXXFLAGS) -o $@ $<

calc_T_P1400: calc_T_P1400.o
	$(CXX) $(CXXFLAGS) -o $@ $<

calc_s0: calc_s0.o calc_distance.o
	$(CXX) $(CXXFLAGS) -o $@ $^

calc_F1400_rcDA: calc_F1400_rcDA.o calc_distance.o
	$(CXX) $(CXXFLAGS) -o $@ $^

calc_F1400_rcDA.o: calc_F1400_rcDA.cc calc_distance.h
	$(CXX) $(CXXFLAGS) -c $<

z_m_simulation: z_m_simulation.o calc_distance.o
	$(CXX) $(CXXFLAGS) -o $@ $^

z_m_simulation.o: z_m_simulation.cc calc_distance.h
	$(CXX) $(CXXFLAGS) -c $<

transfer_sur: transfer_sur.o
	$(CXX) $(CXXFLAGS) $(FIO_LDFLAGS) $(FIO_LIBS) -o $@ $<

calc_nmz: calc_nmz.o
	$(CXX) $(CXXFLAGS) -o $@ $<

calc_nmz.o: calc_nmz.cc adapt_trapezoid.h
	$(CXX) $(CXXFLAGS) -c $<

draw_line: draw_line.o
	$(CXX) $(CXXFLAGS) -o $@ $<

histogram: histogram.o
	$(CXX) $(CXXFLAGS) -o $@ $<


calc_distance.o: calc_distance.cc calc_distance.h
	$(CXX) $(CXXFLAGS) -c $<


calc_NN0: calc_NN0.o
	$(CC) $(CFLAGS) $(GSL_LIBS) -o $@ $<

calc_Fnu_t: calc_Fnu_t.o
	$(CC) $(CFLAGS) $(GSL_LIBS) -o $@ $<


%.o: %.cc
	$(CXX) $(CXXFLAGS) $(FIO_INC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<


clean:
	-rm -f *.o *~ $(TARGETS)

