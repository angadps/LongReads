# @(#)makefile	2-99
# $Revision:$

.SUFFIXES: .o .C .c .cpp

# FLAGS=  -O2
FLAGS= -g -std=c++0x

SRCDIR = .
BINDIR = .
OBJDIR = obj/

HEADER = read.h  snp.h hmm.h

OBJ = read.o  snp.o hmm.o

OBJ3 = program.o

#===============================================================================


program : $(HEADER) $(OBJ) $(OBJ3)
	   g++ $(FLAGS) -o $(BINDIR)/program $(OBJ) $(OBJ3) -lm

#===============================================================================

#===============================================================================

snp.o: snp.h

read.o: read.h

hmm.o: snp.h read.h hmm.h

program.o: hmm.h read.h snp.h

#===============================================================================

#===============================================================================

.cpp.o: $(HEADERI)
	g++ -c  $(FLAGS) $<
.C.o: $(HEADERI)
	g++ -c  $(FLAGS) $<

#checkout:
#	co -l $(REV) $(ALL)

#checkin:
#	ci $(REV)  $(ALL)

clean:
	rm -f *.o core tags

tags: $(SRC)
	ctags *.c

lint:	rm -f lint.out
	lint *.c > lint.out

tar.all: $(SRC) $(DOC) $(HEADERS) makefile
	tar cvf tar.all $(SRC) $(DOC) $(HEADERS) makefile

#===============================================================================
#===============================================================================
