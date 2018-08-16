PREFIX=/usr/local/

libcoscalc.a:coscalc.o
	ar rv $@ $<

coscalc.o:coscalc.cpp coscalc.hpp
	g++ -c $< -O3

clean:
	rm -f *.o

install: libcoscalc.a
	cp coscalc.hpp $(PREFIX)/include/
	cp libcoscalc.a $(PREFIX)/lib/
