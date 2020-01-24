IGSL=$(shell gsl-config --cflags)
LGSL=$(shell gsl-config --libs)
ILHAPDF=$(shell lhapdf-config --cppflags)
LLHAPDF=$(shell lhapdf-config --ldflags)
IALL=$(ILHAPDF) $(IGSL) 
LALL=$(LGSL) -lgfortran $(LLHAPDF)



main: respack.o
	g++ $(IALL) -o pro respack.o  nlo_aux.o $(LALL)
	rm *.o

clean:
	rm *.o

respack.o: nlo_aux.o
	g++ $(IALL) -c respack.cpp -o respack.o

nlo_aux.o:
	gfortran -c nlo_aux.f -o nlo_aux.o
