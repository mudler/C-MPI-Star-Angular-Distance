DISTFILES=$(wildcard *.c *.pl *.pdf) Makefile
DISTNAME=Di_Giacinto
CFLAGS=-Wall -lm

ALL: clean mpi

dist: $(DISTNAME).tar.gz

$(DISTNAME).tar.gz: $(DISTFILES)
	rm -rfv $(DISTNAME)
	mkdir $(DISTNAME)
	cp -rfv dataset $(DISTNAME)/dataset
	cd $(DISTNAME)
	ln $(DISTFILES) $(DISTNAME)
	tar cfz $(DISTNAME).tar.gz $(DISTNAME)
	rm -rfv $(DISTNAME)
mpi: CC=mpicc
mpi: angdist

clean:
	\rm -rfv angdist *.o *~

distclean: clean
	\rm -rfv $(DISTNAME) $(DISTNAME).tar.gz

test: clean
	make mpi

