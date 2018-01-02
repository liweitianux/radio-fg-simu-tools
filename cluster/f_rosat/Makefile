#
# Makefile for f_rosat
#
# Aaron LI
# 2015/04/12
#

CXX= 	g++
CXXFLAGS= 	-Wall


TARGETS= 	calc_Lrosat

all: $(TARGETS)

calc_Lrosat: calc_Lrosat.o f_rosat.o
	$(CXX) $(CXXFLAGS) -o $@ $^

calc_Lrosat.o: calc_Lrosat.cc f_rosat.h
	$(CXX) $(CXXFLAGS) -c $<

f_rosat.o: f_rosat.cc f_rosat.h fconv_table.inc
	$(CXX) $(CXXFLAGS) -c $<

fconv_table.inc: gen_fconv_table.tcl
	xspec - $<

clean:
	-rm -f *.o *~ $(TARGETS)

