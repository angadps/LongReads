# @(#)makefile	2-99
# $Revision:$

.SUFFIXES: .o .C .c .cpp

# FLAGS=  -O2
FLAGS= -g -pg
HAPFLAGS= -Wl,-rpath,/home/singhann/tools/bamtools-master/lib/

SRCDIR = .
BINDIR = .
OBJDIR = obj/

HEADER = read.h  snp.h hmm.h

OBJ = read.o  snp.o hmm.o

OBJ3 = hapmut.o

#===============================================================================

hapmut : $(HEADER) $(OBJ) $(OBJ3)
	g++ -g $(FLAGS) $(HAPFLAGS) -o $(BINDIR)/hapmut $(OBJ) $(OBJ3) -L/home/singhann/tools/bamtools-master/lib/ -lm -lbamtools

#===============================================================================

#===============================================================================

snp.o: snp.h

read.o: read.h

hmm.o: snp.h read.h hmm.h

hapmut.o: hmm.h read.h snp.h

#===============================================================================

#===============================================================================

.cpp.o: $(HEADERI)
	g++ -c -g -I/home/singhann/tools/bamtools-master/include $(FLAGS) $<
.C.o: $(HEADERI)
	g++ -c -g -I/home/singhann/tools/boost_1_54_0 $(FLAGS) $<

#checkout:
#	co -l $(REV) $(ALL)

#checkin:
#	ci $(REV) $(ALL)

clean:
	rm -f *.o core

tags: $(SRC)
	ctags *.c

lint:	rm -f lint.out
	lint *.c > lint.out

tar.all: $(SRC) $(DOC) $(HEADERS) makefile
	tar cvf tar.all $(SRC) $(DOC) $(HEADERS) makefile

#===============================================================================
#===============================================================================
